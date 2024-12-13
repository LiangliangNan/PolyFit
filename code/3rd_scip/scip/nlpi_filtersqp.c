/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    nlpi_filtersqp.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   filterSQP NLP interface
 * @author  Stefan Vigerske
 *
 * @todo scaling
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/time.h>
#endif

/* fallback to non-thread version for windows, because pthread does not exist */
#if defined(_MSC_VER) && defined(SCIP_THREADSAFE)
#undef SCIP_THREADSAFE
#endif

#ifdef SCIP_THREADSAFE
#include <pthread.h>
#endif

#include "scip/nlpi_filtersqp.h"
#include "scip/nlpioracle.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_nlpi.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_solve.h"
#include "scip/pub_misc.h"

#define NLPI_NAME              "filtersqp"                 /**< short concise name of solver */
#define NLPI_DESC              "Sequential Quadratic Programming trust region solver by R. Fletcher and S. Leyffer" /**< description of solver */
#define NLPI_PRIORITY          -1000                       /**< priority of NLP solver */

#define RANDSEED               26051979      /**< initial random seed */
#define MAXPERTURB             0.01          /**< maximal perturbation of bounds in starting point heuristic */
#define MAXNRUNS               3             /**< maximal number of FilterSQP calls per NLP solve (several calls if increasing workspace or decreasing eps) */
#define WORKSPACEGROWTHFACTOR  2             /**< factor by which to increase workspace */
#define MINEPS                 1e-14         /**< minimal FilterSQP epsilon */
#define OPTTOLFACTOR           0.5           /**< factor to apply to optimality tolerance, because FilterSQP do scaling */

/*
 * Data structures
 */

typedef int fint;
typedef double real;
typedef long ftnlen;

struct SCIP_Time
{
   time_t                     sec;           /**< seconds part of time since epoch */
   long                       usec;          /**< micro-seconds part of time */
};
typedef struct SCIP_Time SCIP_TIME;

struct SCIP_NlpiData
{
   SCIP_RANDNUMGEN*            randnumgen;   /**< random number generator, if we have to make up a starting point */
   SCIP_TIME                   starttime;    /**< time at start of solve */
};

struct SCIP_NlpiProblem
{
   SCIP*                       scip;         /**< SCIP data structure */
   SCIP_NLPIORACLE*            oracle;       /**< Oracle-helper to store and evaluate NLP */

   int                         varssize;     /**< length of variables-related arrays, if allocated */
   int                         conssize;     /**< length of constraints-related arrays, if allocated */

   SCIP_Real*                  initguess;    /**< initial values for primal variables, or NULL if not known, size varssize */
   SCIP_Bool                   warmstart;    /**< whether we could warmstart the next solve */

   SCIP_Real*                  primalvalues; /**< primal values of variables in solution, size varssize */
   SCIP_Real*                  consdualvalues;  /**< dual values of constraints in solution, size conssize */
   SCIP_Real*                  varlbdualvalues; /**< dual values of variable lower bounds in solution, size varssize */
   SCIP_Real*                  varubdualvalues; /**< dual values of variable upper bounds in solution, size varssize */

   SCIP_NLPSOLSTAT             solstat;      /**< solution status from last NLP solve */
   SCIP_NLPTERMSTAT            termstat;     /**< termination status from last NLP solve */
   SCIP_Real                   solvetime;    /**< time spend for last NLP solve */
   int                         niterations;  /**< number of iterations for last NLP solve */

   fint                        istat[14];    /**< integer solution statistics from last FilterSQP call */
   real                        rstat[7];     /**< real solution statistics from last FilterSQP call */
   real                        fmin;         /**< lower bound on objective value */
   SCIP_Real                   maxtime;      /**< time limit */

   /* cached FilterSQP data */
   real*                       x;            /**< variable values, size varssize */
   real*                       c;            /**< constraint value, size conssize */
   real*                       lam;          /**< duals, size varssize + conssize */
   real*                       bl;           /**< variable lower bounds and constraint lhs, size varssize + conssize */
   real*                       bu;           /**< variable upper bounds and constraint rhs, size varssize + conssize */
   real*                       s;            /**< scaling factors, size varssize + conssize */
   char*                       cstype;       /**< constraint linearity, size conssize */
   real*                       a;            /**< gradients values, size la[0]-1 */
   fint*                       la;           /**< gradients indices, size lasize */
   int                         lasize;       /**< length of la array */
   fint*                       hessiannz;    /**< nonzero information about Hessian, size hessiannzsize */
   int                         hessiannzsize;/**< length of hessiannz array */
   real*                       ws;           /**< real workspace, size mxwk */
   fint*                       lws;          /**< integer workspace, size mxiwk */
   fint                        mxwk;         /**< size of real workspace */
   fint                        mxiwk;        /**< size of integer workspace */
   real*                       evalbuffer;   /**< buffer to cache evaluation results before passing it to FilterSQP, size evalbufsize */
   int                         evalbufsize;  /**< size of evaluation buffer */
};

/*
 * Local methods
 */

#ifdef FNAME_LCASE_DECOR
# define F77_FUNC(name,NAME) name ## _
#endif
#ifdef FNAME_UCASE_DECOR
# define F77_FUNC(name,NAME) NAME ## _
#endif
#ifdef FNAME_LCASE_NODECOR
# define F77_FUNC(name,NAME) name
#endif
#ifdef FNAME_UCASE_NODECOR
# define F77_FUNC(name,NAME) NAME
#endif

/** FilterSQP main routine.
 *
 * Array a has length nnza, which is the number of nonzeros in the gradient of the objective and the Jacobian.
 * The first entries of a is the objective gradient, next are the gradients of the constraints.
 *
 * Array la has length lamax, which is at least nnza+m+2.
 * la contains the index information of a row-oriented sparse matrix storage. It stores the number of nonzeros, the column indices, and the row starts:
 * la[0] must be set to nnza+1.
 * la[1]..la[nnza] are indices of the variables corresponding to the entries in a (colidx).
 * la[nnza+1]..la[nnza+1+m] contain the index where each row starts in a and la (rowstart).
 */
void F77_FUNC(filtersqp,FILTERSQP)(
  fint*                  n,                  /**< number of variables */
  fint*                  m,                  /**< number of constraints (excluding simple bounds) */
  fint*                  kmax,               /**< maximum size of null-space (at most n) */
  fint*                  maxa,               /**< maximum number of nonzero entries allowed in Jacobian */
  fint*                  maxf,               /**< maximum size of the filter - typically 100 */
  fint*                  mlp,                /**< maximum level parameter for resolving degeneracy in bqpd - typically 100 */
  fint*                  mxwk,               /**< maximum size of real workspace ws */
  fint*                  mxiwk,              /**< maximum size of integer workspace lws */
  fint*                  iprint,             /**< print flag: 0 is quiet, higher is more */
  fint*                  nout,               /**< output channel - 6 for screen */
  fint*                  ifail,              /**< fail flag and warmstart indicator */
  real*                  rho,                /**< initial trust-region radius - default 10 */
  real*                  x,                  /**< starting point and final solution (array of length n) */
  real*                  c,                  /**< final values of general constraints (array of length m) */
  real*                  f,                  /**< final objective value */
  real*                  fmin,               /**< lower bound on objective value (as termination criteria) */
  real*                  bl,                 /**< lower bounds of variables and constraints (array of length n+m) */
  real*                  bu,                 /**< upper bounds of variables and constraints (array of length n+m) */
  real*                  s,                  /**< scaling factors (array of length n+m) */
  real*                  a,                  /**< objective gradient (always dense) and Jacobian (sparse or dense) entries */
  fint*                  la,                 /**< column indices and length of rows of entries in a (if sparse) or leading dimension of a (if dense) */
  real*                  ws,                 /**< real workspace (array of length mxwk) */
  fint*                  lws,                /**< integer workspace (array of length mxiwk) */
  real*                  lam,                /**< Lagrangian multipliers of simple bounds and general constraints (array of length n+m) */
  char*                  cstype,             /**< indicator whether constraint is linear ('L') or nonlinear ('N') (array of length m) */
  real*                  user,               /**< real workspace passed through to user routines */
  fint*                  iuser,              /**< integer workspace passed through to user routines */
  fint*                  maxiter,            /**< iteration limit for SQP solver */
  fint*                  istat,              /**< integer space for solution statistics (array of length 14) */
  real*                  rstat,              /**< real space for solution statitics (array of length 7) */
  ftnlen                 cstype_len          /**< 1 ??? */
  );

void F77_FUNC(objfun,OBJFUN)(real *x, fint *n, real *f, real *user, fint *iuser,
    fint *errflag);

void F77_FUNC(confun,CONFUN)(real *x, fint *n, fint *m, real *c, real *a, fint *la,
    real *user, fint *iuser, fint *errflag);

/* TODO what are the arguments of this function and does it need to be implemented?
 * it's not in the filterSQP manual, but its an undefined symbol in the filterSQP library
void F77_FUNC(objgrad,OBJGRAD)(fint *, fint *, fint *, real *, real *, fint *, fint
    *, real *, fint *, fint *);
*/
void F77_FUNC(objgrad,OBJGRAD)(void);

void F77_FUNC(gradient,GRADIENT)(fint *n, fint *m, fint *mxa, real *x, real *a, fint *la,
    fint *maxa, real *user, fint *iuser, fint *errflag);

void F77_FUNC(hessian,HESSIAN)(real *x, fint *n, fint *m, fint *phase, real *lam,
    real *ws, fint *lws, real *user, fint *iuser,
    fint *l_hess, fint *li_hess, fint *errflag);

/** common block for problemname */
/*lint -esym(754,char_l,pname,*::char_l,*::pname) */
extern struct
{
   fint char_l;
   char pname[10];
} F77_FUNC(cpname,CPNAME);
/*lint -esym(752,cpname_) */

/** common block for Hessian storage set to 0, i.e. NO Hessian */
/*lint -esym(754,*::phr,*::phc) */
extern struct
{
   fint phl, phr, phc;
} F77_FUNC(hessc,HESSC);
/*lint -esym(754,phr,phc) */

/** common block for upper bound on filter */
extern struct
{
   real ubd, tt;
} F77_FUNC(ubdc,UBDC);

/** common block for infinity & epsilon */
extern struct
{
   real infty, eps;
} F77_FUNC(nlp_eps_inf,NLP_EPS_INF);

/** common block for printing from QP solver */
/*lint -esym(754,*::n_bqpd_calls,*::n_bqpd_prfint) */
extern struct
{
   fint n_bqpd_calls, n_bqpd_prfint;
} F77_FUNC(bqpd_count,BQPD_COUNT);
/*lint -esym(752,bqpd_count_) */
/*lint -esym(754,n_bqpd_calls,n_bqpd_prfint) */

/** common for scaling: scale_mode = 0 (none), 1 (variables), 2 (vars+cons) */
/*lint -esym(754,*::phe) */
extern struct
{
   fint scale_mode, phe;
} F77_FUNC(scalec,SCALEC);
/*lint -esym(754,phe) */

#ifdef SCIP_THREADSAFE
static pthread_mutex_t filtersqpmutex = PTHREAD_MUTEX_INITIALIZER; /*lint !e708*/
#endif

static
SCIP_TIME gettime(void)
{
   SCIP_TIME t;
#ifndef _WIN32
   struct timeval tp; /*lint !e86*/
#endif

#ifdef _WIN32
   t.sec = time(NULL);
   t.usec = 0;

#else
   (void) gettimeofday(&tp, NULL);
   t.sec = tp.tv_sec;
   t.usec = tp.tv_usec;
#endif

   return t;
}

/* gives time since starttime (in problem) */
static
SCIP_Real timeelapsed(
   SCIP_NLPIDATA*        nlpidata            /**< NLPI data */
   )
{
   SCIP_TIME now;

   assert(nlpidata != NULL);

   now = gettime();

   /* now should be after startime */
   assert(now.sec >= nlpidata->starttime.sec);
   assert(now.sec > nlpidata->starttime.sec || now.usec >= nlpidata->starttime.usec);

   return (SCIP_Real)(now.sec - nlpidata->starttime.sec) + 1e-6 * (SCIP_Real)(now.usec - nlpidata->starttime.usec);
}

static
SCIP_Bool timelimitreached(
   SCIP_NLPIDATA*        nlpidata,           /**< NLPI data */
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLPI problem */
   )
{
   if( nlpiproblem->maxtime == SCIP_REAL_MAX )  /*lint !e777 */
      return FALSE;

   return timeelapsed(nlpidata) >= nlpiproblem->maxtime;
}

/** Objective function evaluation */ /*lint -e{715} */
void F77_FUNC(objfun,OBJFUN)(
   real*                 x,                  /**< value of current variables (array of length n) */
   fint*                 n,                  /**< number of variables */
   real*                 f,                  /**< buffer to store value of objective function */
   real*                 user,               /**< user real workspace */
   fint*                 iuser,              /**< user integer workspace */
   fint*                 errflag             /**< set to 1 if arithmetic exception occurs, otherwise 0 */
   )
{  /*lint --e{715} */
   SCIP_NLPIPROBLEM* problem;
   real val;

   problem = (SCIP_NLPIPROBLEM*)(void*)iuser;
   assert(problem != NULL);

   if( SCIPisSolveInterrupted(problem->scip) || timelimitreached((SCIP_NLPIDATA*)(void*)user, problem) )
   {
      SCIPdebugMsg(problem->scip, "interrupted or timelimit reached, issuing arithmetic exception in objfun\n");
      *errflag = 1;
      return;
   }

   *errflag = 1;
   if( SCIPnlpiOracleEvalObjectiveValue(problem->scip, problem->oracle, x, &val) == SCIP_OKAY && SCIPisFinite(val) )
   {
      *errflag = 0;
      *f = val;
   }
   else
   {
      SCIPdebugMsg(problem->scip, "arithmetic exception in objfun\n");
   }
}

/** Constraint functions evaluation */ /*lint -e{715} */
void F77_FUNC(confun,CONFUN)(
   real*                 x,                  /**< value of current variables (array of length n) */
   fint*                 n,                  /**< number of variables */
   fint*                 m,                  /**< number of constraints */
   real*                 c,                  /**< buffer to store values of constraints (array of length m) */
   real*                 a,                  /**< Jacobian matrix entries */
   fint*                 la,                 /**< Jacobian index information */
   real*                 user,               /**< user real workspace */
   fint*                 iuser,              /**< user integer workspace */
   fint*                 errflag             /**< set to 1 if arithmetic exception occurs, otherwise 0 */
   )
{  /*lint --e{715} */
   SCIP_NLPIPROBLEM* problem;
   real val;
   int j;

   problem = (SCIP_NLPIPROBLEM*)(void*)iuser;
   assert(problem != NULL);

   *errflag = 0;
   for( j = 0; j < *m; ++j )
   {
      if( SCIPnlpiOracleEvalConstraintValue(problem->scip, problem->oracle, j, x, &val) != SCIP_OKAY || !SCIPisFinite(val) )
      {
         *errflag = 1;
         SCIPdebugMsg(problem->scip, "arithmetic exception in confun for constraint %d\n", j);
         break;
      }
      c[j] = val;
   }
}

/** Objective gradient and Jacobian evaluation
 *
 * \note If an arithmetic exception occurred, then the gradients must not be modified.
 */ /*lint -e{715} */
void
F77_FUNC(gradient,GRADIENT)(
   fint*                 n,                  /**< number of variables */
   fint*                 m,                  /**< number of constraints */
   fint*                 mxa,                /**< actual number of entries in a */
   real*                 x,                  /**< value of current variables (array of length n) */
   real*                 a,                  /**< Jacobian matrix entries */
   fint*                 la,                 /**< Jacobian index information: column indices and pointers to start of each row */
   fint*                 maxa,               /**< maximal size of a */
   real*                 user,               /**< user real workspace */
   fint*                 iuser,              /**< user integer workspace */
   fint*                 errflag             /**< set to 1 if arithmetic exception occurs, otherwise 0 */
   )
{  /*lint --e{715} */
   SCIP_NLPIPROBLEM* problem;
   SCIP_Real dummy;

   problem = (SCIP_NLPIPROBLEM*)(void*)iuser;
   assert(problem != NULL);
   assert(problem->evalbuffer != NULL);
   assert(problem->evalbufsize >= *maxa);

   *errflag = 1;

   if( SCIPnlpiOracleEvalObjectiveGradient(problem->scip, problem->oracle, x, TRUE, &dummy, problem->evalbuffer) == SCIP_OKAY )
   {
      if( SCIPnlpiOracleEvalJacobian(problem->scip, problem->oracle, x, TRUE, NULL, problem->evalbuffer+*n) == SCIP_OKAY )
      {
         BMScopyMemoryArray(a, problem->evalbuffer, *maxa);
         *errflag = 0;
      }
      else
      {
         SCIPdebugMsg(problem->scip, "arithmetic exception in gradient for constraints\n");
      }
   }
   else
   {
      SCIPdebugMsg(problem->scip, "arithmetic exception in gradient for objective\n");
   }
}

/* Objective gradient evaluation */
/*
void F77_FUNC(objgrad,OBJGRAD)(
   fint*,
   fint*,
   fint*,
   real*,
   real*,
   fint*,
   fint*,
   real*,
   fint*,
   fint*
   )
*/
void F77_FUNC(objgrad,OBJGRAD)(void)
{
   SCIPerrorMessage("Objgrad not implemented. Should not be called.\n");
}

/** Hessian of the Lagrangian evaluation
 *
 * phase = 1 : Hessian of the Lagrangian without objective Hessian
 *
 * phase = 2 : Hessian of the Lagrangian (including objective Hessian)
 *
 * \note If an arithmetic exception occurred, then the Hessian must not be modified.
 */ /*lint -e{715} */
void
F77_FUNC(hessian,HESSIAN)(
   real*                 x,                  /**< value of current variables (array of length n) */
   fint*                 n,                  /**< number of variables */
   fint*                 m,                  /**< number of constraints */
   fint*                 phase,              /**< indicates what kind of Hessian matrix is required */
   real*                 lam,                /**< Lagrangian multipliers (array of length n+m) */
   real*                 ws,                 /**< real workspace for Hessian, passed to Wdotd */
   fint*                 lws,                /**< integer workspace for Hessian, passed to Wdotd */
   real*                 user,               /**< user real workspace */
   fint*                 iuser,              /**< user integer workspace */
   fint*                 l_hess,             /**< space of Hessian real storage ws. On entry: maximal space allowed, on exit: actual amount used */
   fint*                 li_hess,            /**< space of Hessian integer storage lws. On entry: maximal space allowed, on exit: actual amount used */
   fint*                 errflag             /**< set to 1 if arithmetic exception occurs, otherwise 0 */
   )
{  /*lint --e{715} */
   SCIP_NLPIPROBLEM* problem;
   SCIP_Real* lambda;
   int nnz;
   int i;

   problem = (SCIP_NLPIPROBLEM*)(void*)iuser;
   assert(problem != NULL);
   assert(problem->evalbuffer != NULL);

   nnz = problem->hessiannz[0]-1;
   assert(problem->evalbufsize >= nnz);

   *errflag = 1;

   /* initialize lambda to -lam */
   BMSallocMemoryArray(&lambda, *m);
   for( i = 0; i < *m; ++i )
      lambda[i] = -lam[*n+i];

   if( SCIPnlpiOracleEvalHessianLag(problem->scip, problem->oracle, x, TRUE, TRUE, (*phase == 1) ? 0.0 : 1.0, lambda, problem->evalbuffer) == SCIP_OKAY )
   {
      *l_hess = nnz;

      BMScopyMemoryArray(ws, problem->evalbuffer, nnz);

      *errflag = 0;

      /* copy the complete problem->hessiannz into lws */
      for( i = 0; i < nnz + *n + 2; ++i )
         lws[i] = problem->hessiannz[i];
      *li_hess = nnz + *n + 2;
   }
   else
   {
      SCIPdebugMsg(problem->scip, "arithmetic exception in hessian\n");
   }

   BMSfreeMemoryArray(&lambda);
}



static
SCIP_RETCODE setupGradients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< NLPI oracle */
   fint**                la,                 /**< buffer to store pointer to sparsity structure */
   int*                  lasize,             /**< buffer to store length of *la array */
   real**                a                   /**< buffer to store pointer to value buffer */
   )
{
   const int* offset;
   const int* col;
   int nnz;  /* number of nonzeros in Jacobian */
   int nvars;
   int ncons;
   int i;
   int c;

   assert(la != NULL);
   assert(lasize != NULL);
   assert(a != NULL);
   assert(*la == NULL);
   assert(*a == NULL);

   nvars = SCIPnlpiOracleGetNVars(oracle);
   ncons = SCIPnlpiOracleGetNConstraints(oracle);

   /* get jacobian sparsity in oracle format: offset are rowstarts in col and col are column indices */
   SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(scip, oracle, &offset, &col) );
   nnz = offset[ncons];

   /* la stores both column indices (first) and rowstarts (at the end) of the objective gradient and Jacobian
    * For the objective gradient, always n entries are taken, the Jacobian is stored sparse
    * la(0) = n+nnz+1 position where rowstarts start in la
    * la(j) column index of objective gradient or Jacobian row, rowwise
    * la(la(0)) position of first nonzero element for objective gradient in a()
    * la(la(0)+i) position of first nonzero element for constraint i gradient in a(), i=1..m
    * la(la(0)+m+1) = n+nnz first unused position in a
    * where n = nvars and m = ncons
    */

   /* need space for la(0) and column indices and rowstarts (1+ncons+1 for objective, constraints, and end (offsets[ncons])) */
   *lasize = 1 + (nvars+nnz) + 1+ncons + 1;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, la, *lasize) );

   (*la)[0] = nvars+nnz+1;

   /* the objective gradient is stored in sparse form */
   for( i = 0; i < nvars; ++i )
      (*la)[1+i] = 1+i;  /* shift by 1 for Fortran */
   (*la)[(*la)[0]] = 1;  /* objective entries start at the beginning in a, shift by 1 for Fortran */

   /* column indicies are as given by col */
   for( i = 0; i < nnz; ++i )
      (*la)[1+nvars+i] = col[i] + 1;  /* shift by 1 for Fortran */

   /* rowstarts are as given by offset, plus extra for objective gradient */
   for( c = 0; c <= ncons; ++c )
      (*la)[(*la)[0]+1+c] = offset[c] + nvars + 1;  /* shift by nvars for objective, shift by 1 for Fortran */

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, a, nvars + nnz) );

#ifdef SCIP_MORE_DEBUG   /* enable for debugging Jacobian */
   for( i = 0; i < 1 + (nvars+nnz) + 1+ncons + 1; ++i )
      printf("la[%2d] = %2d\n", i, (*la)[i]);
#endif

   return SCIP_OKAY;
}

static
SCIP_RETCODE setupHessian(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< NLPI oracle */
   fint**                la,                 /**< buffer to store pointer to Hessian sparsity structure */
   int*                  lasize              /**< buffer to store length of *la array */
   )
{
   const int* offset;
   const int* col;
   int nnz;  /* number of nonzeros in Jacobian */
   int nvars;
   int i;
   int v;

   assert(la != NULL);
   assert(lasize != NULL);
   assert(*la == NULL);

   nvars = SCIPnlpiOracleGetNVars(oracle);

   /* get Hessian sparsity in oracle format: offset are rowstarts in col and col are column indices */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(scip, oracle, &offset, &col) );
   nnz = offset[nvars];

   /* la stores both column indices (first) and rowstarts (at the end) of the (sparse) Hessian
    * la(0) = nnz+1 position where rowstarts start in la
    * la(j) column index of Hessian row, rowwise
    * la(la(0)+i) position of first nonzero element for row i, i = 0..n-1
    * la(la(0)+n) = nnz first unused position in Hessian
    * where n = nvars
    */

   /* need space for la(0) and column indices and rowstarts
    * 1 for la(0)
    * nnz for column indices
    * nvars for rowstarts
    * 1 for first unused position
    */
   *lasize = 1 + nnz + nvars + 1;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, la, *lasize) );

   (*la)[0] = nnz+1;

   /* column indicies are as given by col */
   for( i = 0; i < nnz; ++i )
      (*la)[1+i] = col[i] + 1;  /* shift by 1 for Fortran */

   /* rowstarts are as given by offset */
   for( v = 0; v <= nvars; ++v )
      (*la)[(*la)[0]+v] = offset[v] + 1;  /* shift by 1 for Fortran */

   F77_FUNC(hessc,HESSC).phl = 1;

#ifdef SCIP_MORE_DEBUG /* enable for debugging Hessian */
   for( i = 0; i < 1 + nnz + nvars + 1; ++i )
      printf("lw[%2d] = %2d\n", i, (*la)[i]);
#endif

   return SCIP_OKAY;
}

/** setup starting point for FilterSQP */
static
SCIP_RETCODE setupStart(
   SCIP_NLPIDATA*        data,               /**< NLPI data */
   SCIP_NLPIPROBLEM*     problem,            /**< NLPI problem */
   real*                 x,                  /**< array to store initial values */
   SCIP_Bool*            success             /**< whether we could setup a point in which functions could be evaluated */
   )
{
   int i;
   int n;
   SCIP_Real val;

   assert(data != NULL);
   assert(problem != NULL);
   assert(x != NULL);
   assert(success != NULL);

   n = SCIPnlpiOracleGetNVars(problem->oracle);

   /* setup starting point */
   if( problem->initguess != NULL )
   {
      for( i = 0; i < n; ++i )
         x[i] = problem->initguess[i];
   }
   else
   {
      SCIP_Real lb;
      SCIP_Real ub;

      SCIPdebugMsg(problem->scip, "FilterSQP started without initial primal values; make up something by projecting 0 onto variable bounds and perturb\n");

      if( data->randnumgen == NULL )
      {
         SCIP_CALL( SCIPcreateRandom(problem->scip, &data->randnumgen, RANDSEED, TRUE) );
      }

      for( i = 0; i < n; ++i )
      {
         lb = SCIPnlpiOracleGetVarLbs(problem->oracle)[i];
         ub = SCIPnlpiOracleGetVarUbs(problem->oracle)[i];
         if( lb > 0.0 )
            x[i] = SCIPrandomGetReal(data->randnumgen, lb, lb + MAXPERTURB*MIN(1.0, ub-lb));
         else if( ub < 0.0 )
            x[i] = SCIPrandomGetReal(data->randnumgen, ub - MAXPERTURB*MIN(1.0, ub-lb), ub);
         else
            x[i] = SCIPrandomGetReal(data->randnumgen, MAX(lb, -MAXPERTURB*MIN(1.0, ub-lb)), MIN(ub, MAXPERTURB*MIN(1.0, ub-lb)));
      }
   }

   /* check whether objective and constraints can be evaluated and differentiated once in starting point
    * NOTE: this does not check whether hessian can be computed!
    */
   *success = SCIPnlpiOracleEvalObjectiveValue(problem->scip, problem->oracle, x, &val) == SCIP_OKAY && SCIPisFinite(val);
   *success &= SCIPnlpiOracleEvalObjectiveGradient(problem->scip, problem->oracle, x, FALSE, &val, problem->evalbuffer) == SCIP_OKAY;  /*lint !e514*/
   i = 0;
   for( ; *success && i < SCIPnlpiOracleGetNConstraints(problem->oracle); ++i )
      *success = SCIPnlpiOracleEvalConstraintValue(problem->scip, problem->oracle, i, x, &val) == SCIP_OKAY && SCIPisFinite(val);
   *success &= SCIPnlpiOracleEvalJacobian(problem->scip, problem->oracle, x, FALSE, NULL, problem->evalbuffer) == SCIP_OKAY;  /*lint !e514*/

   if( !*success )
   {
      SCIPdebugMsg(problem->scip, "could not evaluate or constraint %d in %s starting point or Jacobian not available\n", i-1, problem->initguess != NULL ? "provided" : "made up");

      if( problem->initguess != NULL )
      {
         /* forget given starting point and try to make up our own */
         SCIPfreeBlockMemoryArray(problem->scip, &problem->initguess, problem->varssize);
         SCIP_CALL( setupStart(data, problem, x, success) );
      }
   }

   return SCIP_OKAY;
}

/** sets the solstat and termstat to unknown and other, resp. */
static
void invalidateSolution(
   SCIP_NLPIPROBLEM*     problem             /**< data structure of problem */
   )
{
   assert(problem != NULL);

   problem->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->termstat = SCIP_NLPTERMSTAT_OTHER;
}

/** store NLP solve parameters in nlpiproblem */
static
SCIP_RETCODE handleNlpParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     nlpiproblem,        /**< NLP */
   const SCIP_NLPPARAM   param               /**< solve parameters */
   )
{
   assert(scip != NULL);
   assert(nlpiproblem != NULL);

   if( param.fastfail )
   {
      SCIPdebugMsg(scip, "fast fail parameter not supported by FilterSQP interface yet. Ignored.\n");
   }

   nlpiproblem->fmin = param.lobjlimit;

   nlpiproblem->maxtime = param.timelimit;

   return SCIP_OKAY;
}

/** processes results from FilterSQP call */
static
SCIP_RETCODE processSolveOutcome(
   SCIP_NLPIDATA*        nlpidata,           /**< NLPI data */
   SCIP_NLPIPROBLEM*     problem,            /**< NLPI problem */
   fint                  ifail,              /**< fail flag from FilterSQP call */
   SCIP_Real             feastol,            /**< feasibility tolerance */
   SCIP_Real             opttol,             /**< optimality tolerance */
   real*                 x,                  /**< primal solution values from FilterSQP call, or NULL if stopped before filtersqp got called */
   real*                 lam                 /**< dual solution values from FilterSQP call, or NULL if stopped before filtersqp got called */
   )
{
   int i;
   int nvars;
   int ncons;

   assert(problem != NULL);
   assert(ifail >= 0);
   assert((x != NULL) == (lam != NULL));

   problem->solvetime = timeelapsed(nlpidata);

   nvars = SCIPnlpiOracleGetNVars(problem->oracle);
   ncons = SCIPnlpiOracleGetNConstraints(problem->oracle);

   if( ifail <= 8 && x != NULL )
   {
      /* FilterSQP terminated somewhat normally -> store away solution */

      /* make sure we have memory for solution */
      if( problem->primalvalues == NULL )
      {
         assert(problem->varssize >= nvars); /* ensured in nlpiAddVariables */
         assert(problem->conssize >= ncons); /* ensured in nlpiAddConstraints */
         assert(problem->consdualvalues == NULL);  /* if primalvalues == NULL, then also consdualvalues should be NULL */
         assert(problem->varlbdualvalues == NULL); /* if primalvalues == NULL, then also varlbdualvalues should be NULL */
         assert(problem->varubdualvalues == NULL); /* if primalvalues == NULL, then also varubdualvalues should be NULL */

         SCIP_CALL( SCIPallocBlockMemoryArray(problem->scip, &problem->primalvalues, problem->varssize) );
         SCIP_CALL( SCIPallocBlockMemoryArray(problem->scip, &problem->consdualvalues, problem->conssize) );
         SCIP_CALL( SCIPallocBlockMemoryArray(problem->scip, &problem->varlbdualvalues, problem->varssize) );
         SCIP_CALL( SCIPallocBlockMemoryArray(problem->scip, &problem->varubdualvalues, problem->varssize) );
      }
      else
      {
         assert(problem->consdualvalues != NULL);
         assert(problem->varlbdualvalues != NULL);
         assert(problem->varubdualvalues != NULL);
      }

      for( i = 0; i < nvars; ++i )
      {
         problem->primalvalues[i] = x[i];

         /* if dual from FilterSQP is positive, then it belongs to the lower bound, otherwise to the upper bound */
         problem->varlbdualvalues[i] = MAX(0.0,  lam[i]);
         problem->varubdualvalues[i] = MAX(0.0, -lam[i]);
      }

      /* duals from FilterSQP are negated */
      for( i = 0; i < ncons; ++i )
         problem->consdualvalues[i] = -lam[nvars + i];
   }

   /* translate ifail to solution and termination status and decide whether we could warmstart next */
   problem->warmstart = FALSE;
   switch( ifail )
   {
      case 0: /* successful run, solution found */
         assert(problem->rstat[4] <= feastol); /* should be feasible */
         problem->solstat = (problem->rstat[0] <= opttol ? SCIP_NLPSOLSTAT_LOCOPT : SCIP_NLPSOLSTAT_FEASIBLE);
         problem->termstat = SCIP_NLPTERMSTAT_OKAY;
         problem->warmstart = TRUE;
         break;
      case 1: /* unbounded, feasible point with f(x) <= fmin */
         assert(problem->rstat[4] <= feastol); /* should be feasible */
         problem->solstat = SCIP_NLPSOLSTAT_UNBOUNDED;
         if( problem->fmin == SCIP_REAL_MIN )  /*lint !e777*/
            problem->termstat = SCIP_NLPTERMSTAT_OKAY;  /* fmin was not set */
         else
            problem->termstat = SCIP_NLPTERMSTAT_LOBJLIMIT;
         break;
      case 2: /* linear constraints are inconsistent */
         problem->solstat = SCIP_NLPSOLSTAT_GLOBINFEASIBLE;
         problem->termstat =  SCIP_NLPTERMSTAT_OKAY;
         break;
      case 3: /* (locally) nonlinear infeasible, minimal-infeasible solution found */
         /* problem->solstat = (problem->rstat[0] <= opttol ? SCIP_NLPSOLSTAT_LOCINFEASIBLE : SCIP_NLPSOLSTAT_UNKNOWN); */
         problem->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;  /* TODO FilterSQP does not set rstat[0] in this case, assuming local infeasibility is valid */
         problem->termstat =  SCIP_NLPTERMSTAT_OKAY;
         problem->warmstart = TRUE;
        break;
      case 4: /* terminate at point with h(x) <= eps (constraint violation below epsilon) but QP infeasible */
         assert(problem->rstat[4] <= feastol); /* should be feasible */
         problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         problem->termstat =  SCIP_NLPTERMSTAT_NUMERICERROR;
         problem->warmstart = TRUE;
         break;
      case 5: /* termination with rho < eps (trust region radius below epsilon) */
         if( problem->rstat[4] <= feastol )
            problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_NUMERICERROR;
         problem->warmstart = TRUE;
         break;
      case 6: /* termination with iter > max_iter */
         if( problem->rstat[4] <= feastol )
            problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_ITERLIMIT;
         problem->warmstart = TRUE;
         break;
      case 7: /* crash in user routine (IEEE error) could not be resolved, or timelimit reached, or interrupted */
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         if( problem->solvetime >= problem->maxtime )
         {
            problem->termstat =  SCIP_NLPTERMSTAT_TIMELIMIT;
            problem->warmstart = TRUE;
         }
         else if( SCIPisSolveInterrupted(problem->scip) )
            problem->termstat =  SCIP_NLPTERMSTAT_INTERRUPT;
         else
            problem->termstat =  SCIP_NLPTERMSTAT_EVALERROR;
         break;
      case 8: /* unexpect ifail from QP solver */
         if( problem->rstat[4] <= feastol )
            problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_OTHER;
         break;
      case 9: /* not enough REAL workspace */
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_OUTOFMEMORY;
         break;
      case 10: /* not enough INTEGER workspace */
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_OUTOFMEMORY;
         break;
      default:
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_OTHER;
         break;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of NLP solver interface
 */

/** copy method of NLP interface (called when SCIP copies plugins) */
static
SCIP_DECL_NLPICOPY(nlpiCopyFilterSQP)
{
   SCIP_CALL( SCIPincludeNlpSolverFilterSQP(scip) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data */
static
SCIP_DECL_NLPIFREE(nlpiFreeFilterSQP)
{
   assert(nlpi != NULL);
   assert(nlpidata != NULL);
   assert(*nlpidata != NULL);

   if( (*nlpidata)->randnumgen != NULL )
   {
      SCIPfreeRandom(scip, &(*nlpidata)->randnumgen);
   }

   SCIPfreeBlockMemory(scip, nlpidata);
   assert(*nlpidata == NULL);

   return SCIP_OKAY;
}

/** creates a problem instance */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemFilterSQP)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, problem) );
   (*problem)->scip = scip;

   SCIP_CALL( SCIPnlpiOracleCreate(scip, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName(scip, (*problem)->oracle, name) );

   (*problem)->fmin = SCIP_REAL_MIN;
   (*problem)->maxtime = SCIP_REAL_MAX;

   invalidateSolution(*problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** free a problem instance */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemFilterSQP)
{
   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(scip, &(*problem)->oracle) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->initguess, (*problem)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->primalvalues, (*problem)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->consdualvalues, (*problem)->conssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->varlbdualvalues, (*problem)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->varubdualvalues, (*problem)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->cstype, (*problem)->conssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->s, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->bu, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->bl, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->x, (*problem)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->c, (*problem)->conssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->lam, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->a, (*problem)->la != NULL ? (*problem)->la[0]-1 : 0);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->la, (*problem)->lasize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->hessiannz, (*problem)->hessiannzsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->lws, (*problem)->mxiwk);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->ws, (*problem)->mxwk);
   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->evalbuffer, (*problem)->evalbufsize);

   SCIPfreeBlockMemory(scip, problem);
   assert(*problem == NULL);

   return SCIP_OKAY;
}  /*lint !e715*/

/** add variables */
static
SCIP_DECL_NLPIADDVARS(nlpiAddVarsFilterSQP)
{
   int oldnvars;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   oldnvars = SCIPnlpiOracleGetNVars(problem->oracle);

   SCIP_CALL( SCIPnlpiOracleAddVars(scip, problem->oracle, nvars, lbs, ubs, varnames) );

   SCIPfreeBlockMemoryArrayNull(scip, &problem->initguess, problem->varssize);
   invalidateSolution(problem);
   problem->warmstart = FALSE;

   /* increase variables-related arrays in problem, if necessary */
   if( problem->varssize < SCIPnlpiOracleGetNVars(problem->oracle) )
   {
      int newsize = SCIPcalcMemGrowSize(scip, SCIPnlpiOracleGetNVars(problem->oracle));

      if( problem->primalvalues != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->primalvalues, problem->varssize, newsize) );
      }

      if( problem->varlbdualvalues != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->varlbdualvalues, problem->varssize, newsize) );
      }

      if( problem->varubdualvalues != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->varubdualvalues, problem->varssize, newsize) );
      }

      if( problem->x != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->x, problem->varssize, newsize) );
      }

      if( problem->lam != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->lam, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
      }

      if( problem->bl != NULL )
      {
         assert(problem->bu != NULL);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->bl, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->bu, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
      }

      if( problem->s != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->s, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
      }

      problem->varssize = newsize;
   }

   /* update variable bounds in FilterSQP data */
   if( problem->bl != NULL )
   {
      int i;
      int nconss;

      nconss = SCIPnlpiOracleGetNConstraints(problem->oracle);

      /* bl and bu have first variable bounds and then constraint sides
       * copy the constraint sides to their new position before putting in the new variable bounds
       */
      for( i = nconss-1; i >= 0; --i )
      {
         problem->bl[oldnvars+nvars+i] = problem->bl[oldnvars+i];
         problem->bu[oldnvars+nvars+i] = problem->bu[oldnvars+i];
      }

      /* set bounds of new variable */
      for( i = 0; i < nvars; ++i )
      {
         problem->bl[oldnvars+i] = lbs[i];
         problem->bu[oldnvars+i] = ubs[i];
      }
   }

   /* gradients information is out of date now (objective gradient is stored in dense form) */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   SCIPfreeBlockMemoryArrayNull(scip, &problem->la, problem->lasize);

   /* Hessian information is out of date now (no new entries in Hessian, but also empty cols shows up in sparsity info) */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/


/** add constraints */
static
SCIP_DECL_NLPIADDCONSTRAINTS(nlpiAddConstraintsFilterSQP)
{
   int oldnconss;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   oldnconss = SCIPnlpiOracleGetNConstraints(problem->oracle);

   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, problem->oracle, nconss, lhss, rhss, nlininds, lininds, linvals, exprs, names) );

   invalidateSolution(problem);
   problem->warmstart = FALSE;

   /* increase constraints-related arrays in problem, if necessary */
   if( SCIPnlpiOracleGetNConstraints(problem->oracle) > problem->conssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, SCIPnlpiOracleGetNConstraints(problem->oracle));

      if( problem->consdualvalues != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->consdualvalues, problem->conssize, newsize) );
      }

      if( problem->c != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->c, problem->conssize, newsize) );
      }

      if( problem->lam != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->lam, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
      }

      if( problem->bl != NULL )
      {
         assert(problem->bu != NULL);
         assert(problem->cstype != NULL);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->bl, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->bu, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->cstype, problem->conssize, newsize) );
      }

      if( problem->s != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &problem->s, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
      }

      problem->conssize = newsize;
   }

   /* update constraint sides and type in FilterSQP data */
   if( problem->bl != NULL )
   {
      int i;
      int nvars;

      nvars = SCIPnlpiOracleGetNVars(problem->oracle);

      for( i = 0; i < nconss; ++i )
      {
         problem->bl[nvars+oldnconss+i] = lhss[i];
         problem->bu[nvars+oldnconss+i] = rhss[i];
         problem->cstype[oldnconss+i] = SCIPnlpiOracleIsConstraintNonlinear(problem->oracle, oldnconss+i) ? 'N' : 'L';
      }
   }

   /* gradients information is out of date now */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   SCIPfreeBlockMemoryArrayNull(scip, &problem->la, problem->lasize);

   /* Hessian information is out of date now */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected */
static
SCIP_DECL_NLPISETOBJECTIVE( nlpiSetObjectiveFilterSQP )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleSetObjective(scip, problem->oracle,
         constant, nlins, lininds, linvals, expr) );

   invalidateSolution(problem);

   /* gradients info (la,a) should still be ok, as objective gradient is stored in dense form */

   /* Hessian information is out of date now */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** change variable bounds */
static
SCIP_DECL_NLPICHGVARBOUNDS(nlpiChgVarBoundsFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(scip, problem->oracle, nvars, indices, lbs, ubs) );

   invalidateSolution(problem);

   /* update bounds in FilterSQP data */
   if( problem->bl != NULL )
   {
      int i;

      for( i = 0; i < nvars; ++i )
      {
         problem->bl[indices[i]] = lbs[i];
         problem->bu[indices[i]] = ubs[i];
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change constraint bounds */
static
SCIP_DECL_NLPICHGCONSSIDES(nlpiChgConsSidesFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgConsSides(scip, problem->oracle, nconss, indices, lhss, rhss) );

   invalidateSolution(problem);

   /* update constraint sense in FilterSQP data */
   if( problem->bl != NULL )
   {
      int i;
      int nvars;

      nvars = SCIPnlpiOracleGetNVars(problem->oracle);

      for( i = 0; i < nconss; ++i )
      {
         problem->bl[nvars+indices[i]] = lhss[i];
         problem->bu[nvars+indices[i]] = rhss[i];
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** delete a set of variables */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(scip, problem->oracle, dstats) );

   /* @TODO keep initguess and bl, bu for remaining variables? */

   SCIPfreeBlockMemoryArrayNull(scip, &problem->initguess, problem->varssize);
   invalidateSolution(problem);
   problem->warmstart = FALSE;

   SCIPfreeBlockMemoryArrayNull(scip, &problem->bl, problem->varssize + problem->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->bu, problem->varssize + problem->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->cstype, problem->conssize);  /* because we assume that cstype is allocated iff bl is allocated */

   /* gradients information is out of date now (objective gradient is stored in dense form) */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   SCIPfreeBlockMemoryArrayNull(scip, &problem->la, problem->lasize);

   /* Hessian information is out of date now */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** delete a set of constraints */
static
SCIP_DECL_NLPIDELCONSSET(nlpiDelConstraintSetFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(scip, problem->oracle, dstats) );

   invalidateSolution(problem);
   problem->warmstart = FALSE;

   SCIPfreeBlockMemoryArrayNull(scip, &problem->bl, problem->varssize + problem->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->bu, problem->varssize + problem->conssize);  /*lint !e776 */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->cstype, problem->conssize);  /* because we assume that cstype is allocated iff bl is allocated */

   /* gradients information is out of date now */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   SCIPfreeBlockMemoryArrayNull(scip, &problem->la, problem->lasize);

   /* Hessian information is out of date now */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(scip, problem->oracle, idx, nvals, varidxs, vals) );

   invalidateSolution(problem);

   /* gradients information (la,a) may have changed if elements were added or removed
    * (we only care that sparsity doesn't change, not about actual values in a)
    * TODO free only if coefficients were added or removed (SCIPnlpiOracleChgLinearCoefs() could give feedback)
    */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   SCIPfreeBlockMemoryArrayNull(scip, &problem->la, problem->lasize);

   /* Hessian information should still be ok */

   return SCIP_OKAY;
}  /*lint !e715*/

/** replaces the expression of a constraint or objective */
static
SCIP_DECL_NLPICHGEXPR(nlpiChgExprFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExpr(scip, problem->oracle, idxcons, expr) );

   invalidateSolution(problem);

   /* update constraint linearity in FilterSQP data, as we might have changed from linear to nonlinear now */
   if( problem->cstype != NULL && idxcons >= 0 )
      problem->cstype[idxcons] = expr != NULL ? 'N' : 'L';

   /* gradients information (la,a) may have changed */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   SCIPfreeBlockMemoryArrayNull(scip, &problem->la, problem->lasize);

   /* Hessian information may have changed */
   SCIPfreeBlockMemoryArrayNull(scip, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** change the constant offset in the objective */
static
SCIP_DECL_NLPICHGOBJCONSTANT(nlpiChgObjConstantFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(scip, problem->oracle, objconstant) );

   invalidateSolution(problem);

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets initial guess for primal variables */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessFilterSQP)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   if( primalvalues != NULL )
   {
      if( problem->initguess == NULL )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->initguess, problem->varssize) );
      }
      assert(SCIPnlpiOracleGetNVars(problem->oracle) <= problem->varssize);
      BMScopyMemoryArray(problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle));
   }
   else
   {
      SCIPfreeBlockMemoryArrayNull(scip, &problem->initguess, problem->varssize);
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** tries to solve NLP */
static
SCIP_DECL_NLPISOLVE(nlpiSolveFilterSQP)
{
   SCIP_NLPIDATA* data;
   SCIP_Bool success;
   fint n;
   fint m;
   fint kmax;
   fint maxa;
   fint maxf;
   fint mlp;
   fint lh1;
   fint nout;
   fint ifail;
   fint maxiter;
   fint iprint;
   real rho;
   real f;
   real* user;
   fint* iuser;
   ftnlen cstype_len = 1;
   fint minmxwk;
   fint minmxiwk;
   int nruns;
   int i;

   SCIPdebugMsg(scip, "solve with parameters " SCIP_NLPPARAM_PRINT(param));

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   SCIP_CALL( SCIPnlpiOracleResetEvalTime(scip, problem->oracle) );

   if( param.timelimit == 0.0 )
   {
      /* there is nothing we can do if we are not given any time */
      problem->niterations = 0;
      problem->solvetime = 0.0;
      problem->termstat = SCIP_NLPTERMSTAT_TIMELIMIT;
      problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;

      return SCIP_OKAY;
   }

   /* start measuring time */
   data->starttime = gettime();

   SCIP_CALL( handleNlpParam(scip, problem, param) );

   iprint = param.verblevel;

   /* if warmstart parameter is disabled, then we will not warmstart */
   if( !param.warmstart )
      problem->warmstart = FALSE;

   n = SCIPnlpiOracleGetNVars(problem->oracle);
   m = SCIPnlpiOracleGetNConstraints(problem->oracle);
   kmax = n;    /* maximal nullspace dimension */
   maxf = 100;  /* maximal filter length */
   mlp = 100;   /* maximum level of degeneracy */

   /* TODO eventually, the output should be redirected to the message handler,
    * but even to just redirect to some other file, we would have to open the output-unit in Fortran
    */
   nout = 6;   /* output to stdout for now */
   ifail = problem->warmstart ? -1 : 0;  /* -1 for warmstart, otherwise 0 */
   rho = 10.0; /* initial trust-region radius */

   user = (real*)data;
   iuser = (fint*)problem;
   if( problem->warmstart )  /* if warmstart, then need to keep istat[0] */
      memset(problem->istat+1, 0, sizeof(problem->istat)-sizeof(*problem->istat));
   else
      memset(problem->istat, 0, sizeof(problem->istat));
   memset(problem->rstat, 0, sizeof(problem->rstat));
   problem->niterations = 0;

   if( problem->x == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->x, problem->varssize) );
   }
   if( problem->c == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->c, problem->conssize) );
   }
   if( problem->lam == NULL )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &problem->lam, problem->varssize + problem->conssize) );  /*lint !e776 */
   }
   else
   {
      BMSclearMemoryArray(problem->lam, problem->varssize + problem->conssize);
   }
   if( problem->s == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->s, problem->varssize + problem->conssize) );
   }

   if( problem->la == NULL )
   {
      /* allocate la, a and initialize la for Objective Gradient and Jacobian */
      SCIP_CALL( setupGradients(scip, problem->oracle, &problem->la, &problem->lasize, &problem->a) );
   }
   /* maximal number entries in a = nvars+nnz */
   maxa = problem->la[0]-1;

   if( problem->hessiannz == NULL )
   {
      /* allocate and initialize problem->hessiannz for Hessian */
      SCIP_CALL( setupHessian(scip, problem->oracle, &problem->hessiannz, &problem->hessiannzsize) );
   }

   /* setup variable bounds, constraint sides, and constraint types */
   if( problem->bl == NULL )
   {
      assert(problem->bu == NULL);
      assert(problem->cstype == NULL);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->bl, problem->varssize + problem->conssize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->bu, problem->varssize + problem->conssize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->cstype, problem->conssize) );

      BMScopyMemoryArray(problem->bl, SCIPnlpiOracleGetVarLbs(problem->oracle), n);
      BMScopyMemoryArray(problem->bu, SCIPnlpiOracleGetVarUbs(problem->oracle), n);
      for( i = 0; i < m; ++i )
      {
         problem->bl[n+i] = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
         problem->bu[n+i] = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);
         problem->cstype[i] = SCIPnlpiOracleIsConstraintNonlinear(problem->oracle, i) ? 'N' : 'L';
      }
   }

   /* buffer for evaluation results (used in setupStart already) */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &problem->evalbuffer, &problem->evalbufsize, MAX3(n, problem->hessiannz[0], maxa)) );

   /* setup starting point */
   SCIP_CALL( setupStart(data, problem, problem->x, &success) );
   if( !success )
   {
      /* FilterSQP would crash if starting point cannot be evaluated, so give up */
      SCIP_CALL( processSolveOutcome(data, problem, 7, param.feastol, param.opttol, NULL, NULL) );
      return SCIP_OKAY;
   }

   /* setup workspace */
   /* initial guess of real workspace size */
   /* FilterSQP manual: mxwk = 21*n + 8*m + mlp + 8*maxf + kmax*(kmax+9)/2 + nprof, with nprof = 20*n as a good guess */
   /* Bonmin:           mxwk = 21*n + 8*m + mlp + 8*maxf + lh1 + kmax*(kmax+9)/2 + mxwk0,
    *                      with lh1 = nnz_h+8+2*n+m and mxwk0 = 2000000 (parameter) */
   lh1 = problem->hessiannz[0]-1 + 8 + 2*n + m;
   minmxwk = 21*n + 8*m + mlp + 8*maxf + lh1 + kmax*(kmax+9)/2 + MAX(20*n,2000);
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &problem->ws, &problem->mxwk, minmxwk) );

   /* initial guess of integer workspace size */
   /* FilterSQP manual: mxiwk = 13*n + 4*m + mlp + 100 + kmax */
   /* Bonmin:           mxiwk = 13*n + 4*m + mlp + lh1 + kmax + 113 + mxiwk0, with mxiwk0 = 500000 (parameter) */
   minmxiwk = 13*n + 4*m + mlp + lh1 + 100 + kmax + 113 + MAX(5*n,5000);
   if( !problem->warmstart ) /* if warmstart, then lws should remain untouched (n and m didn't change anyway) */
   {
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &problem->lws, &problem->mxiwk, minmxiwk) );
   }
   assert(problem->lws != NULL);
   /* in case of some evalerrors, not clearing ws could lead to valgrind warnings about use of uninitialized memory */
   BMSclearMemoryArray(problem->ws, problem->mxwk);

   /* from here on we are not thread-safe: if intended for multithread use, then protect filtersqp call with mutex
    * NOTE: we need to make sure that we do not return from nlpiSolve before unlocking the mutex
    */
#ifdef SCIP_THREADSAFE
   (void) pthread_mutex_lock(&filtersqpmutex);
#endif

   /* initialize global variables from filtersqp */
   /* FilterSQP eps is tolerance for both feasibility and optimality, and also for trust-region radius, etc. */
   F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps = MIN(param.feastol, param.opttol * OPTTOLFACTOR);
   F77_FUNC(nlp_eps_inf,NLP_EPS_INF).infty = SCIPinfinity(scip);
   F77_FUNC(ubdc,UBDC).ubd = 100.0;
   F77_FUNC(ubdc,UBDC).tt = 1.25;
   F77_FUNC(scalec,SCALEC).scale_mode = 0;

   for( nruns = 1; ; ++nruns )
   {
      maxiter = param.iterlimit - problem->niterations;

      F77_FUNC(filtersqp,FILTERSQP)(
         &n, &m, &kmax, &maxa,
         &maxf, &mlp, &problem->mxwk, &problem->mxiwk,
         &iprint, &nout, &ifail, &rho,
         problem->x, problem->c, &f, &problem->fmin, problem->bl,
         problem->bu, problem->s, problem->a, problem->la, problem->ws,
         problem->lws, problem->lam, problem->cstype, user,
         iuser, &maxiter, problem->istat,
         problem->rstat, cstype_len);

      problem->niterations += problem->istat[1];

      assert(ifail <= 10);
      /* if ifail >= 8 (probably the workspace was too small), then retry with larger workspace
       * if ifail == 0 (local optimal), but absolute violation of KKT too large, then retry with small eps
       */
      if( ifail < 8 && (ifail != 0 || problem->rstat[0] <= param.opttol) )
         break;

      if( param.verblevel > 0 )
      {
         SCIPinfoMessage(scip, NULL, "FilterSQP terminated with status %d in run %d, absolute KKT violation is %g\n", ifail, nruns, problem->rstat[0]);
      }

      /* if iteration or time limit exceeded or solve is interrupted, then don't retry */
      if( problem->niterations >= param.iterlimit || SCIPisSolveInterrupted(scip) || timelimitreached(data, problem) )
      {
         if( param.verblevel > 0 )
         {
            SCIPinfoMessage(scip, NULL, "Time or iteration limit reached or interrupted, not retrying\n");
         }
         break;
      }

      /* if maximal number of runs reached, then stop */
      if( nruns >= MAXNRUNS )
      {
         if( param.verblevel > 0 )
         {
            SCIPinfoMessage(scip, NULL, "Run limit reached, not retrying\n");
         }
         break;
      }

      if( ifail == 0 )
      {
         SCIP_Real epsfactor;

         if( F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps <= MINEPS )
         {
            if( param.verblevel > 0 )
            {
               SCIPinfoMessage(scip, NULL, "Already reached minimal epsilon, not retrying\n");
            }
            break;
         }

         epsfactor = param.opttol / problem->rstat[0];
         assert(epsfactor < 1.0); /* because of the if's above */
         epsfactor *= OPTTOLFACTOR;

         F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps = MAX(MINEPS, F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps * epsfactor);
         if( param.verblevel > 0 )
         {
            SCIPinfoMessage(scip, NULL, "Continue with eps = %g\n", F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps);
         }
         ifail = -1;  /* do warmstart */

         continue;
      }

      /* increase real workspace, if ifail = 9 (real workspace too small) or ifail = 8 (unexpected ifail from QP solver, often also when workspace too small) */
      if( ifail == 8 || ifail == 9 )
      {
         int newsize = SCIPcalcMemGrowSize(scip, WORKSPACEGROWTHFACTOR*problem->mxwk);
         if( BMSreallocBlockMemoryArray(SCIPblkmem(scip), &problem->ws, problem->mxwk, newsize) == NULL )
         {
            /* realloc failed: give up NLP solve */
            problem->mxwk = 0;
            break;
         }
         problem->mxwk = newsize;
      }

      /* increase integer workspace, if ifail = 10 (integer workspace too small) or ifail = 8 (unexpected ifail from QP solver, often also when workspace too small) */
      if( ifail == 8 || ifail == 10 )
      {
         int newsize = SCIPcalcMemGrowSize(scip, WORKSPACEGROWTHFACTOR*problem->mxiwk);
         if( BMSreallocBlockMemoryArray(SCIPblkmem(scip), &problem->lws, problem->mxiwk, newsize) == NULL )
         {
            /* realloc failed: give up NLP solve */
            problem->mxiwk = 0;
            break;
         }
         problem->mxiwk = newsize;

         /* better don't try warmstart for the next trial; warmstart requires that lws is untouched, does extending count as touching? */
         ifail = 0;
      }

      /* reset starting point, in case it was overwritten by failed solve (return can only be SCIP_OKAY, because randnumgen must exist already)
       * NOTE: If initguess is NULL (no user-given starting point), then this will result in a slightly different starting point as in the previous setupStart() call (random numbers)
       *       However, as no point was given, it shouldn't matter which point we actually start from.
       */
      (void) setupStart(data, problem, problem->x, &success);
      assert(success);
   }

#ifdef SCIP_THREADSAFE
   (void) pthread_mutex_unlock(&filtersqpmutex);
#endif

   SCIP_CALL( processSolveOutcome(data, problem, ifail, param.feastol, param.opttol, problem->x, problem->lam) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solution status */
static
SCIP_DECL_NLPIGETSOLSTAT(nlpiGetSolstatFilterSQP)
{
   assert(problem != NULL);

   return problem->solstat;
}  /*lint !e715*/

/** gives termination reason */
static
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetTermstatFilterSQP)
{
   assert(problem != NULL);

   return problem->termstat;
}  /*lint !e715*/

/** gives primal and dual solution values */
static
SCIP_DECL_NLPIGETSOLUTION(nlpiGetSolutionFilterSQP)
{
   assert(problem != NULL);

   if( primalvalues != NULL )
   {
      assert(problem->primalvalues != NULL);

      *primalvalues = problem->primalvalues;
   }

   if( consdualvalues != NULL )
   {
      assert(problem->consdualvalues != NULL);

      *consdualvalues = problem->consdualvalues;
   }

   if( varlbdualvalues != NULL )
   {
      assert(problem->varlbdualvalues != NULL);

      *varlbdualvalues = problem->varlbdualvalues;
   }

   if( varubdualvalues != NULL )
   {
      assert(problem->varubdualvalues != NULL);

      *varubdualvalues = problem->varubdualvalues;
   }

   if( objval != NULL )
   {
      if( problem->primalvalues != NULL )
      {
         /* TODO store last solution value instead of reevaluating the objective function */
         SCIP_CALL( SCIPnlpiOracleEvalObjectiveValue(scip, problem->oracle, problem->primalvalues, objval) );
      }
      else
         *objval = SCIP_INVALID;
   }

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solve statistics */
static
SCIP_DECL_NLPIGETSTATISTICS(nlpiGetStatisticsFilterSQP)
{
   assert(problem != NULL);
   assert(statistics != NULL);

   statistics->niterations = problem->niterations;
   statistics->totaltime = problem->solvetime;
   statistics->evaltime = SCIPnlpiOracleGetEvalTime(scip, problem->oracle);
   statistics->consviol = problem->rstat[4];
   statistics->boundviol = 0.0;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for filterSQP solver and include it into SCIP, if filterSQP is available */
SCIP_RETCODE SCIPincludeNlpSolverFilterSQP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLPIDATA* nlpidata;

   /* create filterSQP solver interface data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &nlpidata) );

   /* create solver interface */
   SCIP_CALL( SCIPincludeNlpi(scip,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyFilterSQP, nlpiFreeFilterSQP, NULL,
         nlpiCreateProblemFilterSQP, nlpiFreeProblemFilterSQP, NULL,
         nlpiAddVarsFilterSQP, nlpiAddConstraintsFilterSQP, nlpiSetObjectiveFilterSQP,
         nlpiChgVarBoundsFilterSQP, nlpiChgConsSidesFilterSQP, nlpiDelVarSetFilterSQP, nlpiDelConstraintSetFilterSQP,
         nlpiChgLinearCoefsFilterSQP, nlpiChgExprFilterSQP,
         nlpiChgObjConstantFilterSQP, nlpiSetInitialGuessFilterSQP, nlpiSolveFilterSQP, nlpiGetSolstatFilterSQP, nlpiGetTermstatFilterSQP,
         nlpiGetSolutionFilterSQP, nlpiGetStatisticsFilterSQP,
         nlpidata) );

   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameFilterSQP(), SCIPgetSolverDescFilterSQP()) );

   return SCIP_OKAY;
}

/** gets string that identifies filterSQP */
const char* SCIPgetSolverNameFilterSQP(
   void
   )
{
   return "filterSQP";  /* TODO version number? */
}

/** gets string that describes filterSQP */
const char* SCIPgetSolverDescFilterSQP(
   void
   )
{
   return NLPI_DESC;
}

/** returns whether filterSQP is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisFilterSQPAvailableFilterSQP(
   void
   )
{
   return TRUE;
}
