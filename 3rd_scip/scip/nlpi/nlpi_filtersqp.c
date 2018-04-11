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

/**@file    nlpi_filtersqp.c
 * @ingroup NLPIS
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

#ifndef NPARASCIP
#include <pthread.h>
#endif

#include "scip/misc.h"
#include "nlpi/nlpi_filtersqp.h"
#include "nlpi/nlpi.h"
#include "nlpi/nlpioracle.h"

#define NLPI_NAME              "filtersqp"                 /* short concise name of solver */
#define NLPI_DESC              "Sequential Quadratic Programming trust region solver by R. Fletcher and S. Leyffer" /* description of solver */
#define NLPI_PRIORITY          -1000                       /* priority of NLP solver */

#define RANDSEED               26051979      /**< initial random seed */
#define MAXPERTURB             0.01          /**< maximal perturbation of bounds in starting point heuristic */
#define MAXNRUNS               3             /**< maximal number of FilterSQP calls per NLP solve (several calls if increasing workspace or decreasing eps) */
#define WORKSPACEGROWTHFACTOR  2             /**< factor by which to increase worksapce */
#define MINEPS                 1e-14         /**< minimal FilterSQP epsilon */
#define OPTTOLFACTOR           0.5           /**< factor to apply to optimality tolerance, because FilterSQP do scaling */
#define DEFAULT_LOBJLIM        (real)(-1e100) /**< default lower objective limit (should mean "unlimited") */
#define DEFAULT_FEASOPTTOL     1e-6          /**< default feasibility and optimality tolerance */
#define DEFAULT_MAXITER        3000          /**< default iteration limit */

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
   BMS_BLKMEM*                 blkmem;       /**< block memory */
   SCIP_MESSAGEHDLR*           messagehdlr;  /**< message handler */
   SCIP_Real                   infinity;     /**< initial value for infinity */
   SCIP_RANDNUMGEN*            randnumgen;   /**< random number generator, if we have to make up a starting point */
   SCIP_TIME                   starttime;    /**< time at start of solve */
};

struct SCIP_NlpiProblem
{
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

   SCIP_Bool                   fromscratch;  /**< value of fromscratch parameter */
   fint                        istat[14];    /**< integer solution statistics from last FilterSQP call */
   real                        rstat[7];     /**< real solution statistics from last FilterSQP call */
   real                        feastol;      /**< user-given feasibility tolerance */
   real                        opttol;       /**< user-given optimality tolerance */
   real                        fmin;         /**< lower bound on objective value */
   SCIP_Real                   maxtime;      /**< time limit */
   fint                        maxiter;      /**< iteration limit */
   fint                        iprint;       /**< print verbosity level */

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
extern struct
{
   fint char_l;
   char pname[10];
} F77_FUNC(cpname,CPNAME);
/*lint -esym(752,cpname_) -esym(754,char_l) -esym(754,pname) */

/** common block for Hessian storage set to 0, i.e. NO Hessian */
extern struct
{
   fint phl, phr, phc;
} F77_FUNC(hessc,HESSC);
/*lint -esym(754,phr) -esym(754,phc) */

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
extern struct
{
   fint n_bqpd_calls, n_bqpd_prfint;
} F77_FUNC(bqpd_count,BQPD_COUNT);
/*lint -esym(752,bqpd_count_) -esym(754,n_bqpd_calls) -esym(754,n_bqpd_prfint) */

/** common for scaling: scale_mode = 0 (none), 1 (variables), 2 (vars+cons) */
extern struct
{
   fint scale_mode, phe;
} F77_FUNC(scalec,SCALEC);
/*lint -esym(754,phe) */

#ifndef NPARASCIP
static pthread_mutex_t filtersqpmutex = PTHREAD_MUTEX_INITIALIZER;
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
   if( nlpiproblem->maxtime == DBL_MAX )  /*lint !e777 */
      return FALSE;

   return timeelapsed(nlpidata) >= nlpiproblem->maxtime;
}

/** Objective function evaluation */
void F77_FUNC(objfun,OBJFUN)(
   real*                 x,                  /**< value of current variables (array of length n) */
   fint*                 n,                  /**< number of variables */
   real*                 f,                  /**< buffer to store value of objective function */
   real*                 user,               /**< user real workspace */
   fint*                 iuser,              /**< user integer workspace */
   fint*                 errflag             /**< set to 1 if arithmetic exception occurs, otherwise 0 */
   )
{ /*lint --e{715} */
   SCIP_NLPIPROBLEM* problem;
   real val;

   problem = (SCIP_NLPIPROBLEM*)(void*)iuser;
   assert(problem != NULL);

   if( timelimitreached((SCIP_NLPIDATA*)(void*)user, problem) )
   {
      SCIPdebugMessage("timelimit reached, issuing arithmetic exception in objfun\n");
      *errflag = 1;
      return;
   }

   *errflag = 1;
   if( SCIPnlpiOracleEvalObjectiveValue(problem->oracle, x, &val) == SCIP_OKAY && SCIPisFinite(val) )
   {
      *errflag = 0;
      *f = val;
   }
   else
   {
      SCIPdebugMessage("arithmetic exception in objfun\n");
   }
}

/** Constraint functions evaluation */
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
{ /*lint --e{715} */
   SCIP_NLPIPROBLEM* problem;
   real val;
   int j;

   problem = (SCIP_NLPIPROBLEM*)(void*)iuser;
   assert(problem != NULL);

   *errflag = 0;
   for( j = 0; j < *m; ++j )
   {
      if( SCIPnlpiOracleEvalConstraintValue(problem->oracle, j, x, &val) != SCIP_OKAY || !SCIPisFinite(val) )
      {
         *errflag = 1;
         SCIPdebugMessage("arithmetic exception in confun for constraint %d\n", j);
         break;
      }
      c[j] = val;
   }
}

/** Objective gradient and Jacobian evaluation
 *
 * \note If an arithmetic exception occurred, then the gradients must not be modified.
 */
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
{ /*lint --e{715} */
   SCIP_NLPIPROBLEM* problem;
   SCIP_Real dummy;

   problem = (SCIP_NLPIPROBLEM*)(void*)iuser;
   assert(problem != NULL);
   assert(problem->evalbuffer != NULL);
   assert(problem->evalbufsize >= *maxa);

   *errflag = 1;

   if( SCIPnlpiOracleEvalObjectiveGradient(problem->oracle, x, TRUE, &dummy, problem->evalbuffer) == SCIP_OKAY )
   {
      if( SCIPnlpiOracleEvalJacobian(problem->oracle, x, TRUE, NULL, problem->evalbuffer+*n) == SCIP_OKAY )
      {
         BMScopyMemoryArray(a, problem->evalbuffer, *maxa);
         *errflag = 0;
      }
      else
      {
         SCIPdebugMessage("arithmetic exception in gradient for constraints\n");
      }
   }
   else
   {
      SCIPdebugMessage("arithmetic exception in gradient for objective\n");
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
 * phase = 2 : Hessian of the Lagrangian (including objective Hessian)
 *
 * \note If an arithmetic exception occurred, then the Hessian must not be modified.
 */
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
{ /*lint --e{715} */
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

   if( SCIPnlpiOracleEvalHessianLag(problem->oracle, x, TRUE, (*phase == 1) ? 0.0 : 1.0, lambda, problem->evalbuffer) == SCIP_OKAY )
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
      SCIPdebugMessage("arithmetic exception in hessian\n");
   }

   BMSfreeMemoryArray(&lambda);
}



static
SCIP_RETCODE setupGradients(
   BMS_BLKMEM*           blkmem,             /**< block memory */
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
   SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(oracle, &offset, &col) );
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
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, la, *lasize) );

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

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, a, nvars + nnz) );

#if 0  /* enable for debugging Jacobian */
   for( i = 0; i < 1 + (nvars+nnz) + 1+ncons + 1; ++i )
      printf("la[%2d] = %2d\n", i, (*la)[i]);
#endif

   return SCIP_OKAY;
}

static
SCIP_RETCODE setupHessian(
   BMS_BLKMEM*           blkmem,             /**< block memory */
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
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(oracle, &offset, &col) );
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
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, la, *lasize) );

   (*la)[0] = nnz+1;

   /* column indicies are as given by col */
   for( i = 0; i < nnz; ++i )
      (*la)[1+i] = col[i] + 1;  /* shift by 1 for Fortran */

   /* rowstarts are as given by offset */
   for( v = 0; v <= nvars; ++v )
      (*la)[(*la)[0]+v] = offset[v] + 1;  /* shift by 1 for Fortran */

   F77_FUNC(hessc,HESSC).phl = 1;

#if 0 /* enable for debugging Hessian */
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

      SCIPdebugMessage("FilterSQP started without initial primal values; make up something by projecting 0 onto variable bounds and perturb\n");

      if( data->randnumgen == NULL )
      {
         SCIP_CALL( SCIPrandomCreate(&data->randnumgen, data->blkmem, RANDSEED) );
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
   *success = SCIPnlpiOracleEvalObjectiveValue(problem->oracle, x, &val) == SCIP_OKAY && SCIPisFinite(val);
   *success &= SCIPnlpiOracleEvalObjectiveGradient(problem->oracle, x, FALSE, &val, problem->evalbuffer) == SCIP_OKAY;  /*lint !e514*/
   i = 0;
   for( ; *success && i < SCIPnlpiOracleGetNConstraints(problem->oracle); ++i )
      *success = SCIPnlpiOracleEvalConstraintValue(problem->oracle, i, x, &val) == SCIP_OKAY && SCIPisFinite(val);
   *success &= SCIPnlpiOracleEvalJacobian(problem->oracle, x, FALSE, NULL, problem->evalbuffer) == SCIP_OKAY;  /*lint !e514*/

   if( !*success )
   {
      SCIPdebugMessage("could not evaluate or constraint %d in %s starting point or Jacobian not available\n", i-1, problem->initguess != NULL ? "provided" : "made up");

      if( problem->initguess != NULL )
      {
         /* forget given starting point and try to make up our own */
         BMSfreeBlockMemoryArray(data->blkmem, &problem->initguess, problem->varssize);
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


/** processes results from FilterSQP call */
static
SCIP_RETCODE processSolveOutcome(
   SCIP_NLPIDATA*        nlpidata,           /**< NLPI data */
   SCIP_NLPIPROBLEM*     problem,            /**< NLPI problem */
   fint                  ifail,              /**< fail flag from FilterSQP call */
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

         SCIP_ALLOC( BMSallocBlockMemoryArray(nlpidata->blkmem, &problem->primalvalues, problem->varssize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(nlpidata->blkmem, &problem->consdualvalues, problem->conssize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(nlpidata->blkmem, &problem->varlbdualvalues, problem->varssize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(nlpidata->blkmem, &problem->varubdualvalues, problem->varssize) );
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
         assert(problem->rstat[4] <= problem->feastol); /* should be feasible */
         problem->solstat = (problem->rstat[0] <= problem->opttol ? SCIP_NLPSOLSTAT_LOCOPT : SCIP_NLPSOLSTAT_FEASIBLE);
         problem->termstat = SCIP_NLPTERMSTAT_OKAY;
         problem->warmstart = TRUE;
         break;
      case 1: /* unbounded, feasible point with f(x) <= fmin */
         assert(problem->rstat[4] <= problem->feastol); /* should be feasible */
         problem->solstat = SCIP_NLPSOLSTAT_UNBOUNDED;
         if( problem->fmin == DEFAULT_LOBJLIM )  /*lint !e777*/
            problem->termstat = SCIP_NLPTERMSTAT_OKAY;  /* fmin was not set */
         else
            problem->termstat = SCIP_NLPTERMSTAT_LOBJLIM;
         break;
      case 2: /* linear constraints are inconsistent */
         problem->solstat = SCIP_NLPSOLSTAT_GLOBINFEASIBLE;
         problem->termstat =  SCIP_NLPTERMSTAT_OKAY;
         break;
      case 3: /* (locally) nonlinear infeasible, minimal-infeasible solution found */
         /* problem->solstat = (problem->rstat[0] <= problem->opttol ? SCIP_NLPSOLSTAT_LOCINFEASIBLE : SCIP_NLPSOLSTAT_UNKNOWN); */
         problem->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;  /* TODO FilterSQP does not set rstat[0] in this case, assuming local infeasibility is valid */
         problem->termstat =  SCIP_NLPTERMSTAT_OKAY;
         problem->warmstart = TRUE;
        break;
      case 4: /* terminate at point with h(x) <= eps (constraint violation below epsilon) but QP infeasible */
         assert(problem->rstat[4] <= problem->feastol); /* should be feasible */
         problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         problem->termstat =  SCIP_NLPTERMSTAT_NUMERR;
         problem->warmstart = TRUE;
         break;
      case 5: /* termination with rho < eps (trust region radius below epsilon) */
         if( problem->rstat[4] <= problem->feastol )
            problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_NUMERR;
         problem->warmstart = TRUE;
         break;
      case 6: /* termination with iter > max_iter */
         if( problem->rstat[4] <= problem->feastol )
            problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_ITLIM;
         problem->warmstart = TRUE;
         break;
      case 7: /* crash in user routine (IEEE error) could not be resolved, or timelimit reached */
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         if( problem->solvetime >= problem->maxtime )
         {
            problem->termstat =  SCIP_NLPTERMSTAT_TILIM;
            problem->warmstart = TRUE;
         }
         else
            problem->termstat =  SCIP_NLPTERMSTAT_EVALERR;
         break;
      case 8: /* unexpect ifail from QP solver */
         if( problem->rstat[4] <= problem->feastol )
            problem->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_OTHER;
         break;
      case 9: /* not enough REAL workspace */
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_MEMERR;
         break;
      case 10: /* not enough INTEGER workspace */
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_MEMERR;
         break;
      default:
         problem->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
         problem->termstat =  SCIP_NLPTERMSTAT_OTHER;
         break;
   }

   return SCIP_OKAY;
}


/** calculate memory size for dynamically allocated arrays (copied from scip/set.c) */
static
int calcGrowSize(
   int                   num                 /**< minimum number of entries to store */
   )
{
   int size;

   /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
   size = 4;
   while( size < num )
      size = (int)(1.2 * size + 4);

   return size;
}


/*
 * Callback methods of NLP solver interface
 */

/** copy method of NLP interface (called when SCIP copies plugins)
 *
 * input:
 *  - blkmem block memory in target SCIP
 *  - sourcenlpi the NLP interface to copy
 *  - targetnlpi buffer to store pointer to copy of NLP interface
 */
static
SCIP_DECL_NLPICOPY( nlpiCopyFilterSQP )
{
   SCIP_NLPIDATA* sourcedata;

   assert(sourcenlpi != NULL);
   assert(targetnlpi != NULL);

   sourcedata = SCIPnlpiGetData(sourcenlpi);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPcreateNlpSolverFilterSQP(blkmem, targetnlpi) );
   assert(*targetnlpi != NULL);

   SCIP_CALL( SCIPnlpiSetRealPar(*targetnlpi, NULL, SCIP_NLPPAR_INFINITY, sourcedata->infinity) );
   SCIP_CALL( SCIPnlpiSetMessageHdlr(*targetnlpi, sourcedata->messagehdlr) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data
 *
 * input:
 *  - nlpi datastructure for solver interface
 */
static
SCIP_DECL_NLPIFREE( nlpiFreeFilterSQP )
{
   SCIP_NLPIDATA* data;

   assert(nlpi != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   if( data->randnumgen != NULL )
   {
      SCIPrandomFree(&data->randnumgen, data->blkmem);
   }

   BMSfreeBlockMemory(data->blkmem, &data);
   assert(data == NULL);

   return SCIP_OKAY;
}

/** gets pointer for NLP solver
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *
 * return: void pointer to solver
 */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerFilterSQP)
{
   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** creates a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer to store the problem data
 *  - name name of problem, can be NULL
 */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemFilterSQP)
{
   SCIP_NLPIDATA* data;

   assert(nlpi    != NULL);
   assert(problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, problem) );
   BMSclearMemory(*problem);

   SCIP_CALL( SCIPnlpiOracleCreate(data->blkmem, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetInfinity((*problem)->oracle, data->infinity) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName((*problem)->oracle, name) );

   (*problem)->feastol = DEFAULT_FEASOPTTOL;
   (*problem)->opttol = DEFAULT_FEASOPTTOL;
   (*problem)->fmin = DEFAULT_LOBJLIM;
   (*problem)->maxtime = DBL_MAX;
   (*problem)->maxiter = DEFAULT_MAXITER;
   (*problem)->iprint = 0;

   invalidateSolution(*problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** free a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer where problem data is stored
 */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemFilterSQP)
{
   SCIP_NLPIDATA* data;

   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(&(*problem)->oracle) );
   }

   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->initguess, (*problem)->varssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->primalvalues, (*problem)->varssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->consdualvalues, (*problem)->conssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->varlbdualvalues, (*problem)->varssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->varubdualvalues, (*problem)->varssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->cstype, (*problem)->conssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->s, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->bu, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->bl, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->x, (*problem)->varssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->c, (*problem)->conssize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->lam, (*problem)->varssize + (*problem)->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->a, (*problem)->la != NULL ? (*problem)->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->la, (*problem)->lasize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->hessiannz, (*problem)->hessiannzsize);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->lws, (*problem)->mxiwk);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->ws, (*problem)->mxwk);
   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->evalbuffer, (*problem)->evalbufsize);

   BMSfreeBlockMemory(data->blkmem, problem);
   assert(*problem == NULL);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets pointer to solver-internal problem instance
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: void pointer to problem instance
 */
static
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerFilterSQP)
{
   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** add variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables
 *  - lbs lower bounds of variables, can be NULL if -infinity
 *  - ubs upper bounds of variables, can be NULL if +infinity
 *  - varnames names of variables, can be NULL
 */
static
SCIP_DECL_NLPIADDVARS( nlpiAddVarsFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;
   int oldnvars;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   oldnvars = SCIPnlpiOracleGetNVars(problem->oracle);

   SCIP_CALL( SCIPnlpiOracleAddVars(problem->oracle, nvars, lbs, ubs, varnames) );

   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->initguess, problem->varssize);
   invalidateSolution(problem);
   problem->warmstart = FALSE;

   /* increase variables-related arrays in problem, if necessary */
   if( problem->varssize < SCIPnlpiOracleGetNVars(problem->oracle) )
   {
      int newsize = calcGrowSize(SCIPnlpiOracleGetNVars(problem->oracle));

      if( problem->primalvalues != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->primalvalues, problem->varssize, newsize) );
      }

      if( problem->varlbdualvalues != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->varlbdualvalues, problem->varssize, newsize) );
      }

      if( problem->varubdualvalues != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->varubdualvalues, problem->varssize, newsize) );
      }

      if( problem->x != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->x, problem->varssize, newsize) );
      }

      if( problem->lam != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->lam, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
      }

      if( problem->bl != NULL )
      {
         assert(problem->bu != NULL);
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->bl, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->bu, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
      }

      if( problem->s != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->s, problem->varssize + problem->conssize, newsize + problem->conssize) );  /*lint !e776*/
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
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian information is out of date now (no new entries in Hessian, but also empty cols shows up in sparsity info) */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/


/** add constraints
 * quadratic coefficients: row oriented matrix for each constraint
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - ncons number of added constraints
 *  - lhss left hand sides of constraints
 *  - rhss right hand sides of constraints
 *  - nlininds number of linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - lininds indices of variables for linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - linvals values of linear coefficient for each constraint
 *    may be NULL in case of no linear part
 *  - nquadrows number of columns in matrix of quadratic part for each constraint
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadrowidxs indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds[.] and quadvals[.]
 *    indices are given w.r.t. quadrowidxs., i.e., quadoffsets[.][i] gives the start index of row quadrowidxs[.][i] in quadvals[.]
 *    quadoffsets[.][nquadrows[.]] gives length of quadinds[.] and quadvals[.]
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadinds column indices w.r.t. quadrowidxs, i.e., quadrowidxs[quadinds[.][i]] gives the index of the variable corresponding
 *    to entry i, entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadvals coefficient values
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    entry of array may be NULL in case of no expression tree
 *    may be NULL in case of no expression tree in any constraint
 *  - exprtrees expression tree for nonquadratic part of constraints
 *    entry of array may be NULL in case of no nonquadratic part
 *    may be NULL in case of no nonquadratic part in any constraint
 *  - names of constraints, may be NULL or entries may be NULL
 */
static
SCIP_DECL_NLPIADDCONSTRAINTS( nlpiAddConstraintsFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;
   int oldnconss;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   oldnconss = SCIPnlpiOracleGetNConstraints(problem->oracle);

   SCIP_CALL( SCIPnlpiOracleAddConstraints(problem->oracle,
         ncons, lhss, rhss,
         nlininds, lininds, linvals,
         nquadelems, quadelems,
         exprvaridxs, exprtrees, names) );

   invalidateSolution(problem);
   problem->warmstart = FALSE;

   /* increase constraints-related arrays in problem, if necessary */
   if( SCIPnlpiOracleGetNConstraints(problem->oracle) > problem->conssize )
   {
      int newsize = calcGrowSize(SCIPnlpiOracleGetNConstraints(problem->oracle));

      if( problem->consdualvalues != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->consdualvalues, problem->conssize, newsize) );
      }

      if( problem->c != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->c, problem->conssize, newsize) );
      }

      if( problem->lam != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->lam, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
      }

      if( problem->bl != NULL )
      {
         assert(problem->bu != NULL);
         assert(problem->cstype != NULL);
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->bl, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->bu, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->cstype, problem->conssize, newsize) );
      }

      if( problem->s != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(nlpidata->blkmem, &problem->s, problem->varssize + problem->conssize, problem->varssize + newsize) );  /*lint !e776*/
      }

      problem->conssize = newsize;
   }

   /* update constraint sides and type in FilterSQP data */
   if( problem->bl != NULL )
   {
      int i;
      int nvars;

      nvars = SCIPnlpiOracleGetNVars(problem->oracle);

      for( i = 0; i < ncons; ++i )
      {
         problem->bl[nvars+oldnconss+i] = lhss[i];
         problem->bu[nvars+oldnconss+i] = rhss[i];
         problem->cstype[oldnconss+i] = SCIPnlpiOracleGetConstraintDegree(problem->oracle, oldnconss+i) <= 1 ? 'L' : 'N';
      }
   }

   /* gradients information is out of date now */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian information is out of date now */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected
 *  May change sparsity pattern.
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nlins number of linear variables
 *  - lininds variable indices
 *    may be NULL in case of no linear part
 *  - linvals coefficient values
 *    may be NULL in case of no linear part
 *  - nquadcols number of columns in matrix of quadratic part
 *  - quadcols indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds and quadvals
 *    quadoffsets[.][nquadcols] gives length of quadinds and quadvals
 *    may be NULL in case of no quadratic part
 *  - quadinds column indices
 *    may be NULL in case of no quadratic part
 *  - quadvals coefficient values
 *    may be NULL in case of no quadratic part
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    may be NULL in case of no expression tree
 *  - exprtree expression tree for nonquadratic part of objective function
 *    may be NULL in case of no nonquadratic part
 *  - constant objective value offset
 */
static
SCIP_DECL_NLPISETOBJECTIVE( nlpiSetObjectiveFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPnlpiOracleSetObjective(problem->oracle,
         constant, nlins, lininds, linvals,
         nquadelems, quadelems,
         exprvaridxs, exprtree) );

   invalidateSolution(problem);

   /* gradients info (la,a) should still be ok, as objective gradient is stored in dense form */

   /* Hessian information is out of date now */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** change variable bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables to change bounds
 *  - indices indices of variables to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
static
SCIP_DECL_NLPICHGVARBOUNDS( nlpiChgVarBoundsFilterSQP )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(problem->oracle, nvars, indices, lbs, ubs) );

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

/** change constraint bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nconss number of constraints to change sides
 *  - indices indices of constraints to change sides
 *  - lhss new left hand sides
 *  - rhss new right hand sides
 */
static
SCIP_DECL_NLPICHGCONSSIDES( nlpiChgConsSidesFilterSQP )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgConsSides(problem->oracle, nconss, indices, lhss, rhss) );

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

/** delete a set of variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 *
 * output:
 *  - dstats new position of var, -1 if var was deleted
 */
static
SCIP_DECL_NLPIDELVARSET( nlpiDelVarSetFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(problem->oracle, dstats) );

   /* @TODO keep initguess and bl, bu for remaining variables? */

   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->initguess, problem->varssize);
   invalidateSolution(problem);
   problem->warmstart = FALSE;

   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->bl, problem->varssize + problem->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->bu, problem->varssize + problem->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->cstype, problem->conssize);  /* because we assume that cstype is allocated iff bl is allocated */

   /* gradients information is out of date now (objective gradient is stored in dense form) */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian information is out of date now */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** delete a set of constraints
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of rows; 1 if row should be deleted, 0 if not
 *
 * output:
 *  - dstats new position of row, -1 if row was deleted
 */
static
SCIP_DECL_NLPIDELCONSSET( nlpiDelConstraintSetFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(problem->oracle, dstats) );

   invalidateSolution(problem);
   problem->warmstart = FALSE;

   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->bl, problem->varssize + problem->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->bu, problem->varssize + problem->conssize);  /*lint !e776 */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->cstype, problem->conssize);  /* because we assume that cstype is allocated iff bl is allocated */

   /* gradients information is out of date now */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian information is out of date now */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint to change
 *  - varidxs indices of variables which coefficient to change
 *  - vals new values for coefficients
 */
static
SCIP_DECL_NLPICHGLINEARCOEFS( nlpiChgLinearCoefsFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(problem->oracle, idx, nvals, varidxs, vals) );

   invalidateSolution(problem);

   /* gradients information (la,a) may have changed if elements were added or removed
    * (we only care that sparsity doesn't change, not about actual values in a)
    * TODO free only if coefficients were added or removed (SCIPnlpiOracleChgLinearCoefs() could give feedback)
    */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian information should still be ok */

   return SCIP_OKAY;
}  /*lint !e715*/

/** changes (or adds) coefficients in the quadratic part of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nentries number of entries in quadratic matrix to change
 *  - rows row indices of entries in quadratic matrix where values should be changed
 *  - cols column indices of entries in quadratic matrix where values should be changed
 *  - values new values for entries in quadratic matrix
 */
static
SCIP_DECL_NLPICHGQUADCOEFS( nlpiChgQuadraticCoefsFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPnlpiOracleChgQuadCoefs(problem->oracle, idx, nquadelems, quadelems) );

   invalidateSolution(problem);

   /* update constraint linearity in FilterSQP data, as we might have changed from linear to nonlinear now */
   if( problem->cstype != NULL && idx >= 0 )
      problem->cstype[idx] = (SCIPnlpiOracleGetConstraintDegree(problem->oracle, idx) <= 1 ? 'L' : 'N');

   /* gradients information (la,a) may have changed if elements were added or removed
    * (we only care that sparsity doesn't change, not about actual values in a)
    * TODO free only if coefficients were added or removed (SCIPnlpiOracleChgLinearCoefs() could give feedback)
    */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian sparsity may have changed if elements were added or removed
    * TODO free only if coefficients were added or removed (SCIPnlpiOracleChgLinearCoefs() could give feedback)
    */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** replaces the expression tree of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp, or NULL
 *  - exprtree new expression tree for constraint or objective, or NULL to only remove previous tree
 */
static
SCIP_DECL_NLPICHGEXPRTREE( nlpiChgExprtreeFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExprtree(problem->oracle, idxcons, exprvaridxs, exprtree) );

   invalidateSolution(problem);

   /* update constraint linearity in FilterSQP data, as we might have changed from linear to nonlinear now */
   if( problem->cstype != NULL && idxcons >= 0 )
      problem->cstype[idxcons] = (SCIPnlpiOracleGetConstraintDegree(problem->oracle, idxcons) <= 1 ? 'L' : 'N');

   /* gradients information (la,a) may have changed */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian information may have changed */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** change one coefficient in the nonlinear part
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - idxparam index of parameter
 *  - value new value for nonlinear parameter
 *
 * return: Error if parameter does not exist
 */
static
SCIP_DECL_NLPICHGNONLINCOEF( nlpiChgNonlinCoefFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExprParam(problem->oracle, idxcons, idxparam, value) );

   invalidateSolution(problem);

   /* gradients information (la,a) may have changed (?) */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->a, problem->la != NULL ? problem->la[0]-1 : 0);
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->la, problem->lasize);

   /* Hessian information may have changed (?) */
   BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->hessiannz, problem->hessiannzsize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** change the constant offset in the objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - objconstant new value for objective constant
 */
static
SCIP_DECL_NLPICHGOBJCONSTANT( nlpiChgObjConstantFilterSQP )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(problem->oracle, objconstant) );

   invalidateSolution(problem);

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets initial guess for primal variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues initial primal values for variables, or NULL to clear previous values
 *  - consdualvalues initial dual values for constraints, or NULL to clear previous values
 *  - varlbdualvalues  initial dual values for variable lower bounds, or NULL to clear previous values
 *  - varubdualvalues  initial dual values for variable upper bounds, or NULL to clear previous values
 */
static
SCIP_DECL_NLPISETINITIALGUESS( nlpiSetInitialGuessFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   if( primalvalues != NULL )
   {
      if( problem->initguess == NULL )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(nlpidata->blkmem, &problem->initguess, problem->varssize) );
      }
      assert(SCIPnlpiOracleGetNVars(problem->oracle) <= problem->varssize);
      BMScopyMemoryArray(problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle));
   }
   else
   {
      BMSfreeBlockMemoryArrayNull(nlpidata->blkmem, &problem->initguess, problem->varssize);
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** tries to solve NLP
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 */
static
SCIP_DECL_NLPISOLVE( nlpiSolveFilterSQP )
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
   real rho;
   real f;
   real* user;
   fint* iuser;
   ftnlen cstype_len = 1;
   fint minmxwk;
   fint minmxiwk;
   int nruns;
   int i;

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   /* start measuring time */
   data->starttime = gettime();

   /* if fromscratch parameter is set, then we will not warmstart */
   if( problem->fromscratch )
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
      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->x, problem->varssize) );
   }
   if( problem->c == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->c, problem->conssize) );
   }
   if( problem->lam == NULL )
   {
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(data->blkmem, &problem->lam, problem->varssize + problem->conssize) );  /*lint !e776 */
   }
   else
   {
      BMSclearMemoryArray(problem->lam, problem->varssize + problem->conssize);
   }
   if( problem->s == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->s, problem->varssize + problem->conssize) );
   }

   if( problem->la == NULL )
   {
      /* allocate la, a and initialize la for Objective Gradient and Jacobian */
      SCIP_CALL( setupGradients(data->blkmem, problem->oracle, &problem->la, &problem->lasize, &problem->a) );
   }
   /* maximal number entries in a = nvars+nnz */
   maxa = problem->la[0]-1;

   if( problem->hessiannz == NULL )
   {
      /* allocate and initialize problem->hessiannz for Hessian */
      SCIP_CALL( setupHessian(data->blkmem, problem->oracle, &problem->hessiannz, &problem->hessiannzsize) );
   }

   /* setup variable bounds, constraint sides, and constraint types */
   if( problem->bl == NULL )
   {
      assert(problem->bu == NULL);
      assert(problem->cstype == NULL);

      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->bl, problem->varssize + problem->conssize) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->bu, problem->varssize + problem->conssize) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->cstype, problem->conssize) );

      BMScopyMemoryArray(problem->bl, SCIPnlpiOracleGetVarLbs(problem->oracle), n);
      BMScopyMemoryArray(problem->bu, SCIPnlpiOracleGetVarUbs(problem->oracle), n);
      for( i = 0; i < m; ++i )
      {
         problem->bl[n+i] = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
         problem->bu[n+i] = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);
         problem->cstype[i] = SCIPnlpiOracleGetConstraintDegree(problem->oracle, i) <= 1 ? 'L' : 'N';
      }
   }

   /* buffer for evaluation results (used in setupStart already) */
   if( problem->evalbufsize < MAX3(n, problem->hessiannz[0], maxa) )
   {
      int newsize = calcGrowSize(MAX3(n, problem->hessiannz[0], maxa));
      SCIP_ALLOC( BMSreallocBlockMemoryArray(data->blkmem, &problem->evalbuffer, problem->evalbufsize, newsize) );
      problem->evalbufsize = newsize;
   }

   /* setup starting point */
   SCIP_CALL( setupStart(data, problem, problem->x, &success) );
   if( !success )
   {
      /* FilterSQP would crash if starting point cannot be evaluated, so give up */
      SCIP_CALL( processSolveOutcome(data, problem, 7, NULL, NULL) );
      return SCIP_OKAY;
   }

   /* setup workspace */
   /* initial guess of real workspace size */
   /* FilterSQP manual: mxwk = 21*n + 8*m + mlp + 8*maxf + kmax*(kmax+9)/2 + nprof, with nprof = 20*n as a good guess */
   /* Bonmin:           mxwk = 21*n + 8*m + mlp + 8*maxf + lh1 + kmax*(kmax+9)/2 + mxwk0,
    *                      with lh1 = nnz_h+8+2*n+m and mxwk0 = 2000000 (parameter) */
   lh1 = problem->hessiannz[0]-1 + 8 + 2*n + m;
   minmxwk = 21*n + 8*m + mlp + 8*maxf + lh1 + kmax*(kmax+9)/2 + MAX(20*n,2000);
   if( problem->ws == NULL )
   {
      problem->mxwk = minmxwk;
      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->ws, problem->mxwk) );
   }
   else if( problem->mxwk < minmxwk )
   {
      int newsize = calcGrowSize(minmxwk);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(data->blkmem, &problem->ws, problem->mxwk, newsize) );
      problem->mxwk = newsize;
   }

   /* initial guess of integer workspace size */
   /* FilterSQP manual: mxiwk = 13*n + 4*m + mlp + 100 + kmax */
   /* Bonmin:           mxiwk = 13*n + 4*m + mlp + lh1 + kmax + 113 + mxiwk0, with mxiwk0 = 500000 (parameter) */
   minmxiwk = 13*n + 4*m + mlp + lh1 + 100 + kmax + 113 + MAX(5*n,5000);
   if( problem->lws == NULL )
   {
      assert(!problem->warmstart);

      problem->mxiwk = minmxiwk;
      SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &problem->lws, problem->mxiwk) );
   }
   else if( problem->mxiwk < minmxiwk && !problem->warmstart ) /* if warmstart, then lws should remain untouched (n and m didn't change anyway) */
   {
      int newsize = calcGrowSize(minmxiwk);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(data->blkmem, &problem->lws, problem->mxiwk, newsize) );
      problem->mxiwk = newsize;
   }
   /* in case of some evalerrors, not clearing ws could lead to valgrind warnings about use of uninitialized memory */
   memset(problem->ws, 0, problem->mxwk * sizeof(real));

   /* from here on we are not thread-safe: if intended for multithread use, then protect filtersqp call with mutex
    * NOTE: we need to make sure that we do not return from nlpiSolve before unlocking the mutex
    */
#ifndef NPARASCIP
   pthread_mutex_lock(&filtersqpmutex);
#endif

   /* initialize global variables from filtersqp */
   /* FilterSQP eps is tolerance for both feasibility and optimality, and also for trust-region radius, etc. */
   F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps = MIN(problem->feastol, problem->opttol * OPTTOLFACTOR);
   F77_FUNC(nlp_eps_inf,NLP_EPS_INF).infty = SCIPnlpiOracleGetInfinity(problem->oracle);
   F77_FUNC(ubdc,UBDC).ubd = 100.0;
   F77_FUNC(ubdc,UBDC).tt = 1.25;
   F77_FUNC(scalec,SCALEC).scale_mode = 0;

   for( nruns = 1; ; ++nruns )
   {
      maxiter = problem->maxiter - problem->niterations;

      F77_FUNC(filtersqp,FILTERSQP)(
         &n, &m, &kmax, &maxa,
         &maxf, &mlp, &problem->mxwk, &problem->mxiwk,
         &problem->iprint, &nout, &ifail, &rho,
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
      if( ifail < 8 && (ifail != 0 || problem->rstat[0] <= problem->opttol) )
         break;

      if( problem->iprint > 0 )
      {
         SCIPmessagePrintInfo(data->messagehdlr, "FilterSQP terminated with status %d in run %d, absolute KKT violation is %g\n", ifail, nruns, problem->rstat[0]);
      }

      /* if iteration or time limit exceeded, then don't retry */
      if( problem->niterations >= problem->maxiter || timelimitreached(data, problem) )
      {
         if( problem->iprint > 0 )
         {
            SCIPmessagePrintInfo(data->messagehdlr, "Time or iteration limit reached, not retrying\n");
         }
         break;
      }

      /* if maximal number of runs reached, then stop */
      if( nruns >= MAXNRUNS )
      {
         if( problem->iprint > 0 )
         {
            SCIPmessagePrintInfo(data->messagehdlr, "Run limit reached, not retrying\n");
         }
         break;
      }

      if( ifail == 0 )
      {
         SCIP_Real epsfactor;

         if( F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps <= MINEPS )
         {
            if( problem->iprint > 0 )
            {
               SCIPmessagePrintInfo(data->messagehdlr, "Already reached minimal epsilon, not retrying\n");
            }
            break;
         }

         epsfactor = problem->opttol / problem->rstat[0];
         assert(epsfactor < 1.0); /* because of the if's above */
         epsfactor *= OPTTOLFACTOR;

         F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps = MAX(MINEPS, F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps * epsfactor);
         if( problem->iprint > 0 )
         {
            SCIPmessagePrintInfo(data->messagehdlr, "Continue with eps = %g\n", F77_FUNC(nlp_eps_inf,NLP_EPS_INF).eps);
         }
         ifail = -1;  /* do warmstart */

         continue;
      }

      /* increase real workspace, if ifail = 9 (real workspace too small) or ifail = 8 (unexpected ifail from QP solver, often also when workspace too small) */
      if( ifail == 8 || ifail == 9 )
      {
         int newsize = calcGrowSize(WORKSPACEGROWTHFACTOR*problem->mxwk);
         if( BMSreallocBlockMemoryArray(data->blkmem, &problem->ws, problem->mxwk, newsize) == NULL )
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
         int newsize = calcGrowSize(WORKSPACEGROWTHFACTOR*problem->mxiwk);
         if( BMSreallocBlockMemoryArray(data->blkmem, &problem->lws, problem->mxiwk, newsize) == NULL )
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

#ifndef NPARASCIP
   pthread_mutex_unlock(&filtersqpmutex);
#endif

   SCIP_CALL( processSolveOutcome(data, problem, ifail, problem->x, problem->lam) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solution status
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Solution Status
 */
static
SCIP_DECL_NLPIGETSOLSTAT( nlpiGetSolstatFilterSQP )
{
   assert(problem != NULL);

   return problem->solstat;
}  /*lint !e715*/

/** gives termination reason
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Termination Status
 */
static
SCIP_DECL_NLPIGETTERMSTAT( nlpiGetTermstatFilterSQP )
{
   assert(problem != NULL);

   return problem->termstat;
}  /*lint !e715*/

/** gives primal and dual solution values
 *
 * solver can return NULL in dual values if not available
 * but if solver provides dual values for one side of variable bounds, then it must also provide those for the other side
 *
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues buffer to store pointer to array to primal values, or NULL if not needed
 *  - consdualvalues buffer to store pointer to array to dual values of constraints, or NULL if not needed
 *  - varlbdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - varubdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - objval buffer store the objective value, or NULL if not needed
 */
static
SCIP_DECL_NLPIGETSOLUTION( nlpiGetSolutionFilterSQP )
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
         SCIP_CALL( SCIPnlpiOracleEvalObjectiveValue(problem->oracle, problem->primalvalues, objval) );
      }
      else
         *objval = SCIP_INVALID;
   }

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives solve statistics
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - statistics pointer to store statistics
 *
 * output:
 *  - statistics solve statistics
 */
static
SCIP_DECL_NLPIGETSTATISTICS( nlpiGetStatisticsFilterSQP )
{
   assert(problem != NULL);

   SCIPnlpStatisticsSetNIterations(statistics, problem->niterations);
   SCIPnlpStatisticsSetTotalTime(statistics, problem->solvetime);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gives required size of a buffer to store a warmstart object
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - size pointer to store required size for warmstart buffer
 *
 * output:
 *  - size required size for warmstart buffer
 */
static
SCIP_DECL_NLPIGETWARMSTARTSIZE( nlpiGetWarmstartSizeFilterSQP )
{
   SCIPerrorMessage("method of filtersqp nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** stores warmstart information in buffer
 *
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer memory to store warmstart information
 *
 * output:
 *  - buffer warmstart information in solver specific data structure
 */
static
SCIP_DECL_NLPIGETWARMSTARTMEMO( nlpiGetWarmstartMemoFilterSQP )
{
   SCIPerrorMessage("method of filtersqp nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets warmstart information in solver
 *
 * write warmstart to buffer
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer warmstart information
 */
static
SCIP_DECL_NLPISETWARMSTARTMEMO( nlpiSetWarmstartMemoFilterSQP )
{
   SCIPerrorMessage("method of filtersqp nonlinear solver is not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival pointer to store the parameter value
 *
 * output:
 *  - ival parameter value
 */
static
SCIP_DECL_NLPIGETINTPAR( nlpiGetIntParFilterSQP )
{
   assert(nlpi != NULL);
   assert(ival != NULL);
   assert(problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      *ival = problem->fromscratch ? 1 : 0;
      break;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      *ival = problem->iprint;
      break;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      SCIPerrorMessage("relative objective tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      SCIPerrorMessage("objective limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      SCIPerrorMessage("infinity parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      *ival = problem->maxiter;
      break;
   }

   case SCIP_NLPPAR_TILIM:
   {
      SCIPerrorMessage("time limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("optfile parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      *ival = 0;
      break;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival parameter value
 */
static
SCIP_DECL_NLPISETINTPAR( nlpiSetIntParFilterSQP )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      if( ival == 0 || ival == 1 )
      {
         problem->fromscratch = (SCIP_Bool)ival;
      }
      else
      {
         SCIPerrorMessage("Value %d for parameter from scratch out of range {0, 1}\n", ival);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      if( ival >= 0 )
      {
         problem->iprint = ival;
      }
      else
      {
         SCIPerrorMessage("Value %d for parameter from verbosity level out of range\n", ival);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      SCIPerrorMessage("relative objective tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      SCIPerrorMessage("objective limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      SCIPerrorMessage("infinity parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      if( ival >= 0 )
      {
         problem->maxiter = ival;
      }
      else
      {
         SCIPerrorMessage("Value %d for parameter iteration limit is negative\n", ival);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_TILIM:
   {
      SCIPerrorMessage("time limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("optfile parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      if( ival == 0 || ival == 1 )
      {
         SCIPdebugMessage("fast fail parameter not supported by FilterSQP interface yet. Ignored.\n");
      }
      else
      {
         SCIPerrorMessage("Value %d for parameter fastfail out of range {0, 1}\n", ival);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to FilterSQP interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval pointer to store the parameter value
 *
 * output:
 *  - dval parameter value
 */
static
SCIP_DECL_NLPIGETREALPAR( nlpiGetRealParFilterSQP )
{
   assert(nlpi != NULL);
   assert(dval != NULL);
   assert(type == SCIP_NLPPAR_INFINITY || problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      SCIPerrorMessage("from scratch parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      SCIPerrorMessage("verbosity level parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      *dval = problem->feastol;
      break;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      *dval = problem->opttol;
      break;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      *dval = problem->fmin;
      break;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      if( problem )
      {
         *dval = SCIPnlpiOracleGetInfinity(problem->oracle);
      }
      else
      {
         SCIP_NLPIDATA* data = SCIPnlpiGetData(nlpi);
         assert(data != NULL);
         *dval = data->infinity;
      }
      break;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      SCIPerrorMessage("iteration limit parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_TILIM:
   {
      *dval = problem->maxtime;
      break;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("option file parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to FilterSQP interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval parameter value
 */
static
SCIP_DECL_NLPISETREALPAR( nlpiSetRealParFilterSQP )
{
   assert(nlpi != NULL);
   assert(type == SCIP_NLPPAR_INFINITY || problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      SCIPerrorMessage("from scratch parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      SCIPerrorMessage("verbosity level parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      if( dval >= 0 )
      {
         problem->feastol = dval;
      }
      else
      {
         SCIPerrorMessage("Value %g for parameter feasibility tolerance is negative\n", dval);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      if( dval >= 0 )
      {
         problem->opttol = dval;
      }
      else
      {
         SCIPerrorMessage("Value %g for parameter relative objective tolerance is negative\n", dval);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      problem->fmin = dval;
      break;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      if( dval < 0.0 )
         return SCIP_PARAMETERWRONGVAL;
      if( problem )
      {
         SCIP_CALL( SCIPnlpiOracleSetInfinity(problem->oracle, dval) );
      }
      else
      {
         SCIP_NLPIDATA* data = SCIPnlpiGetData(nlpi);
         assert(data != NULL);
         data->infinity = dval;
      }
      break;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      SCIPerrorMessage("iteration limit parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_TILIM:
   {
      if( dval >= 0 )
      {
         problem->maxtime = dval;
      }
      else
      {
         SCIPerrorMessage("Value %g for parameter time limit is negative\n", dval);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("option file parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to FilterSQP interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval pointer to store the string value, the user must not modify the string
 *
 * output:
 *  - sval parameter value
 */
static
SCIP_DECL_NLPIGETSTRINGPAR( nlpiGetStringParFilterSQP )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      SCIPerrorMessage("from scratch parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      SCIPerrorMessage("verbosity level parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      SCIPerrorMessage("objective tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      SCIPerrorMessage("objective limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      SCIPerrorMessage("infinity parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      SCIPerrorMessage("iteration limit parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_TILIM:
   {
      SCIPerrorMessage("time limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      *sval = NULL;
      return SCIP_OKAY;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to FilterSQP interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }
}  /*lint !e715*/

/** sets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval parameter value
 */
static
SCIP_DECL_NLPISETSTRINGPAR( nlpiSetStringParFilterSQP )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      SCIPerrorMessage("from scratch parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      SCIPerrorMessage("verbosity level parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      SCIPerrorMessage("objective tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      SCIPerrorMessage("objective limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      SCIPerrorMessage("infinity parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      SCIPerrorMessage("iteration limit parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_TILIM:
   {
      SCIPerrorMessage("time limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPmessagePrintWarning(SCIPnlpiGetData(nlpi)->messagehdlr, "Parameter optfile not supported by FilterSQP interface. Ignored.\n");
      return SCIP_OKAY;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      SCIPerrorMessage("fastfail parameter is of type int.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to Ipopt interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }
}  /*lint !e715*/

/** sets message handler for message output
 *
 * input:
 *  - nlpi NLP interface structure
 *  - messagehdlr SCIP message handler, or NULL to suppress all output
 */
static
SCIP_DECL_NLPISETMESSAGEHDLR( nlpiSetMessageHdlrFilterSQP )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   nlpidata->messagehdlr = messagehdlr;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for FilterSQP solver */
SCIP_RETCODE SCIPcreateNlpSolverFilterSQP(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
   )
{
   SCIP_NLPIDATA* nlpidata;

   assert(blkmem != NULL);
   assert(nlpi   != NULL);

   /* create filterSQP solver interface data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &nlpidata) );

   nlpidata->blkmem = blkmem;
   nlpidata->messagehdlr = NULL;
   nlpidata->infinity = SCIP_DEFAULT_INFINITY;
   nlpidata->randnumgen = NULL;

   /* create solver interface */
   SCIP_CALL( SCIPnlpiCreate(nlpi,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyFilterSQP, nlpiFreeFilterSQP, nlpiGetSolverPointerFilterSQP,
         nlpiCreateProblemFilterSQP, nlpiFreeProblemFilterSQP, nlpiGetProblemPointerFilterSQP,
         nlpiAddVarsFilterSQP, nlpiAddConstraintsFilterSQP, nlpiSetObjectiveFilterSQP,
         nlpiChgVarBoundsFilterSQP, nlpiChgConsSidesFilterSQP, nlpiDelVarSetFilterSQP, nlpiDelConstraintSetFilterSQP,
         nlpiChgLinearCoefsFilterSQP, nlpiChgQuadraticCoefsFilterSQP, nlpiChgExprtreeFilterSQP, nlpiChgNonlinCoefFilterSQP,
         nlpiChgObjConstantFilterSQP, nlpiSetInitialGuessFilterSQP, nlpiSolveFilterSQP, nlpiGetSolstatFilterSQP, nlpiGetTermstatFilterSQP,
         nlpiGetSolutionFilterSQP, nlpiGetStatisticsFilterSQP,
         nlpiGetWarmstartSizeFilterSQP, nlpiGetWarmstartMemoFilterSQP, nlpiSetWarmstartMemoFilterSQP,
         nlpiGetIntParFilterSQP, nlpiSetIntParFilterSQP, nlpiGetRealParFilterSQP, nlpiSetRealParFilterSQP, nlpiGetStringParFilterSQP, nlpiSetStringParFilterSQP,
         nlpiSetMessageHdlrFilterSQP,
         nlpidata) );

   return SCIP_OKAY;
}

/** gets string that identifies filterSQP (version number) */
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
