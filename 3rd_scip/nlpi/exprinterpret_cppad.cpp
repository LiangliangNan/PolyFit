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

/**@file    exprinterpret_cppad.cpp
 * @brief   methods to interpret (evaluate) an expression tree "fast" using CppAD
 * @ingroup EXPRINTS
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "nlpi/pub_expr.h"
#include "nlpi/exprinterpret.h"

#include <cmath>
#include <vector>
using std::vector;

/* Turn off lint warning "747: Significant prototype coercion" and "732: Loss of sign".
 * The first warning is generated for expressions like t[0], where t is a vector, since 0 is an integer constant, but a
 * size_t is expected (usually long unsigned). The second is generated for expressions like t[n], where n is an
 * integer. Both code pieces are likely to be correct. It seems to be impossible to inhibit these messages for
 * vector<*>::operator[] only. */
/*lint --e{747,732}*/

/* Turn off lint info "1702 operator '...' is both an ordinary function 'CppAD::operator...' and a member function 'CppAD::SCIPInterval::operator...'.
 * However, the functions have different signatures (the CppAD working on double, the SCIPInterval member
 * function working on SCIPInterval's.
 */
/*lint --e{1702}*/

/* defining NO_CPPAD_USER_ATOMIC disables the use of our own implementation of derivaties of power operators
 * via CppAD's user-atomic function feature
 * our customized implementation should give better results (tighter intervals) for the interval data type
 */
/* #define NO_CPPAD_USER_ATOMIC */

/** sign of a value (-1 or +1)
 * 
 * 0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

/* in order to use intervals as operands in CppAD,
 * we need to include the intervalarithext.h very early and require the interval operations to be in the CppAD namespace */
#define SCIPInterval_NAMESPACE CppAD
#include "nlpi/intervalarithext.h"

namespace CppAD
{
   SCIP_Real SCIPInterval::infinity = SCIP_DEFAULT_INFINITY;
}
using CppAD::SCIPInterval;

/* CppAD needs to know a fixed upper bound on the number of threads at compile time.
 * It is wise to set it to a power of 2, so that if the tape id overflows, it is likely to start at 0 again, which avoids difficult to debug errors.
 */
#ifndef CPPAD_MAX_NUM_THREADS
#ifndef NPARASCIP
#define CPPAD_MAX_NUM_THREADS 64
#else
#define CPPAD_MAX_NUM_THREADS 1
#endif
#endif

#include <cppad/cppad.hpp>
#include <cppad/utility/error_handler.hpp>

/* CppAD is not thread-safe by itself, but uses some static datastructures
 * To run it in a multithreading environment, a special CppAD memory allocator that is aware of the multiple threads has to be used.
 * This allocator requires to know the number of threads and a thread number for each thread.
 * To implement this, we follow the team_pthread example of CppAD, which uses pthread's thread-specific data management.
 */
#ifndef NPARASCIP
#include <pthread.h>

/** mutex for locking in pthread case */
static pthread_mutex_t cppadmutex = PTHREAD_MUTEX_INITIALIZER;

/** key for accessing thread specific information */
static pthread_key_t thread_specific_key;

/** currently registered number of threads */
static size_t ncurthreads = 0;

/** CppAD callback function that indicates whether we are running in parallel mode */
static
bool in_parallel(void)
{
   return ncurthreads > 1;
}

/** CppAD callback function that returns the number of the current thread
 *
 * assigns a new number to the thread if new
 */
static
size_t thread_num(void)
{
   size_t threadnum;
   void* specific;

   specific = pthread_getspecific(thread_specific_key);

   /* if no data for this thread yet, then assign a new thread number to the current thread
    * we store the thread number incremented by one, to distinguish the absence of data (=0) from existing data
    */
   if( specific == NULL )
   {
      pthread_mutex_lock(&cppadmutex);

      SCIPdebugMessage("Assigning thread number %lu to thread %p.\n", (long unsigned int)ncurthreads, (void*)pthread_self());

      pthread_setspecific(thread_specific_key, (void*)(ncurthreads + 1));

      threadnum = ncurthreads;

      ++ncurthreads;

      pthread_mutex_unlock(&cppadmutex);

      assert(pthread_getspecific(thread_specific_key) != NULL);
      assert((size_t)pthread_getspecific(thread_specific_key) == threadnum + 1);
   }
   else
   {
      threadnum = (size_t)(specific) - 1;
   }

   assert(threadnum < ncurthreads);

   return threadnum;
}

/** sets up CppAD's datastructures for running in multithreading mode
 *
 *  It must be called once before multithreading is started.
 */
static
char init_parallel(void)
{
   pthread_key_create(&thread_specific_key, NULL);

   CppAD::thread_alloc::parallel_setup(CPPAD_MAX_NUM_THREADS, in_parallel, thread_num);
   CppAD::parallel_ad<double>();
   CppAD::parallel_ad<SCIPInterval>();

   return 0;
}

/** a dummy variable that is initialized to the result of init_parallel
 *
 *  The purpose is to make sure that init_parallel() is called before any multithreading is started.
 */
#if !defined(_MSC_VER)
__attribute__ ((unused))
#endif
static char init_parallel_return = init_parallel();

#endif // NPARASCIP

/** definition of CondExpOp for SCIPInterval (required by CppAD) */
inline
SCIPInterval CondExpOp(
   enum CppAD::CompareOp cop,
   const SCIPInterval&   left,
   const SCIPInterval&   right,
   const SCIPInterval&   trueCase,
   const SCIPInterval&   falseCase)
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "SCIPInterval CondExpOp(...)",
      "Error: cannot use CondExp with an interval type"
      );

   return SCIPInterval();
}

/** another function required by CppAD */
inline
bool IdenticalPar(
   const SCIPInterval&   x                   /**< operand */
   )
{  /*lint --e{715}*/
   return true;
}

/** returns whether the interval equals [0,0] */
inline
bool IdenticalZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return (x == 0.0);
}

/** returns whether the interval equals [1,1] */
inline
bool IdenticalOne(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return (x == 1.0);
}

/** yet another function that checks whether two intervals are equal */
inline
bool IdenticalEqualPar(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   return (x == y);
}

/** greater than zero not defined for intervals */
inline
bool GreaterThanZero(
   const SCIPInterval&   x                   /**< operand */
   )
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "GreaterThanZero(x)",
      "Error: cannot use GreaterThanZero with interval"
      );

   return false;
}

/** greater than or equal zero not defined for intervals */
inline
bool GreaterThanOrZero(
   const SCIPInterval&   x                   /**< operand */
   )
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__ ,
      "GreaterThanOrZero(x)",
      "Error: cannot use GreaterThanOrZero with interval"
      );

   return false;
}

/** less than not defined for intervals */
inline
bool LessThanZero(
   const SCIPInterval&   x                   /**< operand */
   )
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "LessThanZero(x)",
      "Error: cannot use LessThanZero with interval"
      );

   return false;
}

/** less than or equal not defined for intervals */
inline
bool LessThanOrZero(
   const SCIPInterval&   x                   /**< operand */
   )
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "LessThanOrZero(x)",
      "Error: cannot use LessThanOrZero with interval"
      );

   return false;
}

/** conversion to integers not defined for intervals */
inline
int Integer(
   const SCIPInterval&   x                   /**< operand */
   )
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "Integer(x)",
      "Error: cannot use Integer with interval"
      );

   return 0;
}

/** absolute zero multiplication
 *
 * @return [0,0] if first argument is [0,0] independent of whether the second argument is an empty interval or not
 */
inline
SCIPInterval azmul(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   if( x.inf == 0.0 && x.sup == 0.0 )
      return SCIPInterval(0.0, 0.0);
   return x * y;
}

/** printing of an interval (required by CppAD) */
inline
std::ostream& operator<<(std::ostream& out, const SCIP_INTERVAL& x)
{
   out << '[' << x.inf << ',' << x.sup << ']';
   return out;
}

using CppAD::AD;

/** expression interpreter */
struct SCIP_ExprInt
{
   BMS_BLKMEM*           blkmem;             /**< block memory data structure */
};

/** expression specific interpreter data */
struct SCIP_ExprIntData
{
public:
   /** constructor */
   SCIP_ExprIntData()
      : val(0.0), need_retape(true), int_need_retape(true), need_retape_always(false), userevalcapability(SCIP_EXPRINTCAPABILITY_ALL), blkmem(NULL), root(NULL)
   { }

   /** destructor */
   ~SCIP_ExprIntData()
   { }/*lint --e{1540}*/

   vector< AD<double> >  X;                  /**< vector of dependent variables */
   vector< AD<double> >  Y;                  /**< result vector */ 
   CppAD::ADFun<double>  f;                  /**< the function to evaluate as CppAD object */

   vector<double>        x;                  /**< current values of dependent variables */
   double                val;                /**< current function value */
   bool                  need_retape;        /**< will retaping be required for the next point evaluation? */

   vector< AD<SCIPInterval> > int_X;         /**< interval vector of dependent variables */
   vector< AD<SCIPInterval> > int_Y;         /**< interval result vector */
   CppAD::ADFun<SCIPInterval> int_f;         /**< the function to evaluate on intervals as CppAD object */

   vector<SCIPInterval>  int_x;              /**< current interval values of dependent variables */
   SCIPInterval          int_val;            /**< current interval function value */
   bool                  int_need_retape;    /**< will retaping be required for the next interval evaluation? */

   bool                  need_retape_always; /**< will retaping be always required? */
   SCIP_EXPRINTCAPABILITY userevalcapability; /**< (intersection of) capabilities of evaluation rountines of user expressions */

   BMS_BLKMEM*           blkmem;             /**< block memory used to allocate expresstion tree */
   SCIP_EXPR*            root;               /**< copy of expression tree; @todo we should not need to make a copy */
};

#ifndef NO_CPPAD_USER_ATOMIC

/** computes sparsity of jacobian for a univariate function during a forward sweep
 *
 *  For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
 *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
 */
static
bool univariate_for_sparse_jac(
   size_t                     q,             /**< number of columns in R */
   const CppAD::vector<bool>& r,             /**< sparsity of R, columnwise */
   CppAD::vector<bool>&       s              /**< vector to store sparsity of S, columnwise */
   )
{
   assert(r.size() == q);
   assert(s.size() == q);

   s = r;

   return true;
}

/** Computes sparsity of jacobian during a reverse sweep
 *
 *  For a q x 1 matrix R, we have to return the sparsity pattern of the q x 1 matrix S(x) = R * f'(x).
 *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
 */
static
bool univariate_rev_sparse_jac(
   size_t                     q,             /**< number of rows in R */
   const CppAD::vector<bool>& r,             /**< sparsity of R, rowwise */
   CppAD::vector<bool>&       s              /**< vector to store sparsity of S, rowwise */
   )
{
   assert(r.size() == q);
   assert(s.size() == q);

   s = r;

   return true;
}

/** computes sparsity of hessian during a reverse sweep
 *
 *  Assume V(x) = (g(f(x)))'' R  with f(x) = x^p for a function g:R->R and a matrix R.
 *  we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
 */
static
bool univariate_rev_sparse_hes(
   const CppAD::vector<bool>& vx,            /**< indicates whether argument is a variable, or empty vector */
   const CppAD::vector<bool>& s,             /**< sparsity pattern of S = g'(y) */
   CppAD::vector<bool>&  t,                  /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
   size_t                q,                  /**< number of columns in R, U, and V */
   const CppAD::vector<bool>& r,             /**< sparsity pattern of R */
   const CppAD::vector<bool>& u,             /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
   CppAD::vector<bool>&  v                   /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
   )
{  /*lint --e{439,715}*/  /* @todo take vx into account */
   assert(r.size() == q);
   assert(s.size() == 1);
   assert(t.size() == 1);
   assert(u.size() == q);
   assert(v.size() == q);

   // T(x) = g'(f(x)) * f'(x) = S * f'(x), and f' is not identically 0
   t[0] = s[0];

   // V(x) = g''(f(x)) f'(x) f'(x) R + g'(f(x)) f''(x) R466
   //      = f'(x) U + S f''(x) R, with f'(x) and f''(x) not identically 0
   v = u;
   if( s[0] )
      for( size_t j = 0; j < q; ++j )
         if( r[j] )
            v[j] = true;

   return true;
}


/** Automatic differentiation of x -> x^p, p>=2 integer, as CppAD user-atomic function.
 *
 *  This class implements forward and reverse operations for the function x -> x^p for use within CppAD.
 *  While CppAD would implement integer powers as a recursion of multiplications, we still use pow functions as they allow us to avoid overestimation in interval arithmetics.
 *
 *  @todo treat the exponent as a (variable) argument to the function, with the assumption that we never differentiate w.r.t. it (this should make the approach threadsafe again)
 */
template<class Type>
class atomic_posintpower : public CppAD::atomic_base<Type>
{
public:
   atomic_posintpower()
   : CppAD::atomic_base<Type>("posintpower"),
     exponent(0)
   {
      /* indicate that we want to use bool-based sparsity pattern */
      this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
   }

private:
   /** exponent value for next call to forward or reverse */
   int exponent;

   /** stores exponent value corresponding to next call to forward or reverse
    *
    * how is this supposed to be threadsafe? (we use only one global instantiation of this class)
    */
   virtual void set_id(size_t id)
   {
      exponent = (int) id;
   }

   /** forward sweep of positive integer power
    *
    * Given the taylor coefficients for x, we have to compute the taylor coefficients for f(x),
    * that is, given tx = (x, x', x'', ...), we compute the coefficients ty = (y, y', y'', ...)
    * in the taylor expansion of f(x) = x^p.
    * Thus, y   = x^p
    *           = tx[0]^p,
    *       y'  = p * x^(p-1) * x'
    *           = p * tx[0]^(p-1) * tx[1],
    *       y'' = 1/2 * p * (p-1) * x^(p-2) * x'^2 + p * x^(p-1) * x''
    *           = 1/2 * p * (p-1) * tx[0]^(p-2) * tx[1]^2 + p * tx[0]^(p-1) * tx[2]
    */
   bool forward(
      size_t                     q,          /**< lowest order Taylor coefficient that we are evaluating */
      size_t                     p,          /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<bool>& vx,         /**< indicates whether argument is a variable, or empty vector */
      CppAD::vector<bool>&       vy,         /**< vector to store which function values depend on variables, or empty vector */
      const CppAD::vector<Type>& tx,         /**< values for taylor coefficients of x */
      CppAD::vector<Type>&       ty          /**< vector to store taylor coefficients of y */
      )
   {
      assert(exponent > 1);
      assert(tx.size() >= p+1);
      assert(ty.size() >= p+1);
      assert(q <= p);

      if( vx.size() > 0 )
      {
         assert(vx.size() == 1);
         assert(vy.size() == 1);
         assert(p == 0);

         vy[0] = vx[0];
      }

      if( q == 0 /* q <= 0 && 0 <= p */ )
      {
         ty[0] = CppAD::pow(tx[0], exponent);
      }

      if( q <= 1 && 1 <= p )
      {
         ty[1] = CppAD::pow(tx[0], exponent-1) * tx[1];
         ty[1] *= double(exponent);
      }

      if( q <= 2 && 2 <= p )
      {
         if( exponent > 2 )
         {
            // ty[2] = 1/2 * exponent * (exponent-1) * pow(tx[0], exponent-2) * tx[1] * tx[1] + exponent * pow(tx[0], exponent-1) * tx[2];
            ty[2]  = CppAD::pow(tx[0], exponent-2) * tx[1] * tx[1];
            ty[2] *= (exponent-1) / 2.0;
            ty[2] += CppAD::pow(tx[0], exponent-1) * tx[2];
            ty[2] *= exponent;
         }
         else
         {
            assert(exponent == 2);
            // ty[2] = 1/2 * exponent * tx[1] * tx[1] + exponent * tx[0] * tx[2];
            ty[2]  = tx[1] * tx[1] + 2.0 * tx[0] * tx[2];
         }
      }

      /* higher order derivatives not implemented */
      if( p > 2 )
         return false;

      return true;
   }

   /** reverse sweep of positive integer power
    *
    * Assume y(x) is a function of the taylor coefficients of f(x) = x^p for x, i.e.,
    *   y(x) = [ x^p, p * x^(p-1) * x', p * (p-1) * x^(p-2) * x'^2 + p * x^(p-1) * x'', ... ].
    * Then in the reverse sweep we have to compute the elements of \f$\partial h / \partial x^[l], l = 0, ..., k,\f$
    * where x^[l] is the l'th taylor coefficient (x, x', x'', ...) and h(x) = g(y(x)) for some function g:R^k -> R.
    * That is, we have to compute
    *\f$
    * px[l] = \partial h / \partial x^[l] = (\partial g / \partial y) * (\partial y / \partial x^[l])
    *       = \sum_{i=0}^k (\partial g / \partial y_i) * (\partial y_i / \partial x^[l])
    *       = \sum_{i=0}^k py[i] * (\partial y_i / \partial x^[l])
    * \f$
    *
    * For k = 0, this means
    *\f$
    * px[0] = py[0] * (\partial y_0 / \partial x^[0])
    *       = py[0] * (\partial x^p / \partial x)
    *       = py[0] * p * tx[0]^(p-1)
    *\f$
    *
    * For k = 1, this means
    * \f$
    * px[0] = py[0] * (\partial y_0 / \partial x^[0]) + py[1] * (\partial y_1 / \partial x^[0])
    *       = py[0] * (\partial x^p / \partial x)     + py[1] * (\partial (p * x^(p-1) * x') / \partial x)
    *       = py[0] * p * tx[0]^(p-1)                 + py[1] * p * (p-1) * tx[0]^(p-2) * tx[1]
    * px[1] = py[0] * (\partial y_0 / \partial x^[1]) + py[1] * (\partial y_1 / \partial x^[1])
    *       = py[0] * (\partial x^p / \partial x')    + py[1] * (\partial (p * x^(p-1) x') / \partial x')
    *       = py[0] * 0                               + py[1] * p * tx[0]^(p-1)
    * \f$
    */
   bool reverse(
      size_t                     p,          /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<Type>& tx,         /**< values for taylor coefficients of x */
      const CppAD::vector<Type>& ty,         /**< values for taylor coefficients of y */
      CppAD::vector<Type>&       px,         /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
      const CppAD::vector<Type>& py          /**< values for partial derivatives of g(x) w.r.t. y */
      )
   { /*lint --e{715}*/
      assert(exponent > 1);
      assert(px.size() >= p+1);
      assert(py.size() >= p+1);
      assert(tx.size() >= p+1);

      switch( p )
      {
         case 0:
            // px[0] = py[0] * exponent * pow(tx[0], exponent-1);
            px[0]  = py[0] * CppAD::pow(tx[0], exponent-1);
            px[0] *= exponent;
            break;

         case 1:
            // px[0] = py[0] * exponent * pow(tx[0], exponent-1) + py[1] * exponent * (exponent-1) * pow(tx[0], exponent-2) * tx[1];
            px[0]  = py[1] * tx[1] * CppAD::pow(tx[0], exponent-2);
            px[0] *= exponent-1;
            px[0] += py[0] * CppAD::pow(tx[0], exponent-1);
            px[0] *= exponent;
            // px[1] = py[1] * exponent * pow(tx[0], exponent-1);
            px[1]  = py[1] * CppAD::pow(tx[0], exponent-1);
            px[1] *= exponent;
            break;

         default:
            return false;
      }

      return true;
   }

   using CppAD::atomic_base<Type>::for_sparse_jac;

   /** computes sparsity of jacobian during a forward sweep
    *
    * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
    * Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool for_sparse_jac(
      size_t                     q,          /**< number of columns in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, columnwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, columnwise */
      )
   {
      return univariate_for_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_jac;

   /** computes sparsity of jacobian during a reverse sweep
    *
    *  For a q x 1 matrix R, we have to return the sparsity pattern of the q x 1 matrix S(x) = R * f'(x).
    *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool rev_sparse_jac(
      size_t                     q,          /**< number of rows in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, rowwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, rowwise */
      )
   {
      return univariate_rev_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_hes;

   /** computes sparsity of hessian during a reverse sweep
    *
    *  Assume V(x) = (g(f(x)))'' R  with f(x) = x^p for a function g:R->R and a matrix R.
    *  we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
    */
   bool rev_sparse_hes(
      const CppAD::vector<bool>&   vx,       /**< indicates whether argument is a variable, or empty vector */
      const CppAD::vector<bool>&   s,        /**< sparsity pattern of S = g'(y) */
      CppAD::vector<bool>&         t,        /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
      size_t                       q,        /**< number of columns in R, U, and V */
      const CppAD::vector<bool>&   r,        /**< sparsity pattern of R */
      const CppAD::vector<bool>&   u,        /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
      CppAD::vector<bool>&         v         /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
      )
   {
      return univariate_rev_sparse_hes(vx, s, t, q, r, u, v);
   }
};

/** power function with natural exponents */
template<class Type>
static
void posintpower(
   const vector<Type>&   in,                 /**< vector which first argument is base */
   vector<Type>&         out,                /**< vector where to store result in first argument */
   size_t                exponent            /**< exponent */
   )
{
   static atomic_posintpower<typename Type::value_type> pip;
   pip(in, out, exponent);
}

#else

/** power function with natural exponents */
template<class Type>
void posintpower(
   const vector<Type>&   in,                 /**< vector which first argument is base */
   vector<Type>&         out,                /**< vector where to store result in first argument */
   size_t                exponent            /**< exponent */
   )
{
   out[0] = pow(in[0], (int)exponent);
}

#endif


#ifndef NO_CPPAD_USER_ATOMIC

/** Automatic differentiation of x -> sign(x)abs(x)^p, p>=1, as CppAD user-atomic function.
 *
 *  This class implements forward and reverse operations for the function x -> sign(x)abs(x)^p for use within CppAD.
 *  While we otherwise would have to use discontinuous sign and abs functions, our own implementation allows to provide
 *  a continuously differentiable function.
 *
 *  @todo treat the exponent as a (variable) argument to the function, with the assumption that we never differentiate w.r.t. it (this should make the approach threadsafe again)
 */
template<class Type>
class atomic_signpower : public CppAD::atomic_base<Type>
{
public:
   atomic_signpower()
   : CppAD::atomic_base<Type>("signpower"),
     exponent(0.0)
   {
      /* indicate that we want to use bool-based sparsity pattern */
      this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
   }

private:
   /** exponent for use in next call to forward or reverse */
   SCIP_Real exponent;

   /** stores exponent corresponding to next call to forward or reverse
    *
    *  How is this supposed to be threadsafe? (we use only one global instantiation of this class)
    */
   virtual void set_id(size_t id)
   {
      exponent = SCIPexprGetSignPowerExponent((SCIP_EXPR*)(void*)id);
   }

   /** forward sweep of signpower
    *
    * Given the taylor coefficients for x, we have to compute the taylor coefficients for f(x),
    * that is, given tx = (x, x', x'', ...), we compute the coefficients ty = (y, y', y'', ...)
    * in the taylor expansion of f(x) = sign(x)abs(x)^p.
    * Thus, y   = sign(x)abs(x)^p
    *           = sign(tx[0])abs(tx[0])^p,
    *       y'  = p * abs(x)^(p-1) * x'
    *           = p * abs(tx[0])^(p-1) * tx[1],
    *       y'' = 1/2 * p * (p-1) * sign(x) * abs(x)^(p-2) * x'^2 + p * abs(x)^(p-1) * x''
    *           = 1/2 * p * (p-1) * sign(tx[0]) * abs(tx[0])^(p-2) * tx[1]^2 + p * abs(tx[0])^(p-1) * tx[2]
    */
   bool forward(
      size_t                      q,         /**< lowest order Taylor coefficient that we are evaluating */
      size_t                      p,         /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<bool>&  vx,        /**< indicates whether argument is a variable, or empty vector */
      CppAD::vector<bool>&        vy,        /**< vector to store which function values depend on variables, or empty vector */
      const CppAD::vector<Type>&  tx,        /**< values for taylor coefficients of x */
      CppAD::vector<Type>&        ty         /**< vector to store taylor coefficients of y */
      )
   {
      assert(exponent > 0.0);
      assert(tx.size() >= p+1);
      assert(ty.size() >= p+1);
      assert(q <= p);

      if( vx.size() > 0 )
      {
         assert(vx.size() == 1);
         assert(vy.size() == 1);
         assert(p == 0);

         vy[0] = vx[0];
      }

      if( q == 0 /* q <= 0 && 0 <= p */ )
      {
         ty[0] = SIGN(tx[0]) * pow(REALABS(tx[0]), exponent);
      }

      if( q <= 1 && 1 <= p )
      {
            ty[1] = pow(REALABS(tx[0]), exponent - 1.0) * tx[1];
            ty[1] *= exponent;
      }

      if( q <= 2 && 2 <= p )
      {
         if( exponent != 2.0 )
         {
            ty[2]  = SIGN(tx[0]) * pow(REALABS(tx[0]), exponent - 2.0) * tx[1] * tx[1];
            ty[2] *= (exponent - 1.0) / 2.0;
            ty[2] += pow(REALABS(tx[0]), exponent - 1.0) * tx[2];
            ty[2] *= exponent;
         }
         else
         {
            // y'' = 2 (1/2 * sign(x) * x'^2 + |x|*x'') = sign(tx[0]) * tx[1]^2 + 2 * abs(tx[0]) * tx[2]
            ty[2]  = SIGN(tx[0]) * tx[1] * tx[1];
            ty[2] += 2.0 * REALABS(tx[0]) * tx[2];
         }
      }

      /* higher order derivatives not implemented */
      if( p > 2 )
         return false;

      return true;
   }

   /** reverse sweep of signpower
    *
    * Assume y(x) is a function of the taylor coefficients of f(x) = sign(x)|x|^p for x, i.e.,
    *   y(x) = [ f(x), f'(x), f''(x), ... ].
    * Then in the reverse sweep we have to compute the elements of \f$\partial h / \partial x^[l], l = 0, ..., k,\f$
    * where x^[l] is the l'th taylor coefficient (x, x', x'', ...) and h(x) = g(y(x)) for some function g:R^k -> R.
    * That is, we have to compute
    *\f$
    * px[l] = \partial h / \partial x^[l] = (\partial g / \partial y) * (\partial y / \partial x^[l])
    *       = \sum_{i=0}^k (\partial g / \partial y_i) * (\partial y_i / \partial x^[l])
    *       = \sum_{i=0}^k py[i] * (\partial y_i / \partial x^[l])
    *\f$
    *
    * For k = 0, this means
    *\f$
    * px[0] = py[0] * (\partial y_0 / \partial x^[0])
    *       = py[0] * (\partial f(x) / \partial x)
    *       = py[0] * p * abs(tx[0])^(p-1)
    * \f$
    *
    * For k = 1, this means
    *\f$
    * px[0] = py[0] * (\partial y_0  / \partial x^[0]) + py[1] * (\partial y_1   / \partial x^[0])
    *       = py[0] * (\partial f(x) / \partial x)     + py[1] * (\partial f'(x) / \partial x)
    *       = py[0] * p * abs(tx[0])^(p-1)             + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
    * px[1] = py[0] * (\partial y_0  / \partial x^[1]) + py[1] * (\partial y_1 / \partial x^[1])
    *       = py[0] * (\partial f(x) / \partial x')    + py[1] * (\partial f'(x) / \partial x')
    *       = py[0] * 0                                + py[1] * p * abs(tx[0])^(p-1)
    * \f$
    */
   bool reverse(
      size_t                      p,         /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<Type>&  tx,        /**< values for taylor coefficients of x */
      const CppAD::vector<Type>&  ty,        /**< values for taylor coefficients of y */
      CppAD::vector<Type>&        px,        /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
      const CppAD::vector<Type>&  py         /**< values for partial derivatives of g(x) w.r.t. y */
      )
   { /*lint --e{715}*/
      assert(exponent > 1);
      assert(px.size() >= p+1);
      assert(py.size() >= p+1);
      assert(tx.size() >= p+1);

      switch( p )
      {
      case 0:
         // px[0] = py[0] * p * pow(abs(tx[0]), p-1);
         px[0]  = py[0] * pow(REALABS(tx[0]), exponent - 1.0);
         px[0] *= p;
         break;

      case 1:
         if( exponent != 2.0 )
         {
            // px[0] = py[0] * p * abs(tx[0])^(p-1) + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
            px[0]  = py[1] * tx[1] * pow(REALABS(tx[0]), exponent - 2.0) * SIGN(tx[0]);
            px[0] *= exponent - 1.0;
            px[0] += py[0] * pow(REALABS(tx[0]), exponent - 1.0);
            px[0] *= exponent;
            // px[1] = py[1] * p * abs(tx[0])^(p-1)
            px[1]  = py[1] * pow(REALABS(tx[0]), exponent - 1.0);
            px[1] *= exponent;
         }
         else
         {
            // px[0] = py[0] * 2.0 * abs(tx[0]) + py[1] * 2.0 * sign(tx[0]) * tx[1]
            px[0]  = py[1] * tx[1] * SIGN(tx[0]);
            px[0] += py[0] * REALABS(tx[0]);
            px[0] *= 2.0;
            // px[1] = py[1] * 2.0 * abs(tx[0])
            px[1]  = py[1] * REALABS(tx[0]);
            px[1] *= 2.0;
         }
         break;

      default:
         return false;
      }

      return true;
   }

   using CppAD::atomic_base<Type>::for_sparse_jac;

   /** computes sparsity of jacobian during a forward sweep
    *
    * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
    * Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool for_sparse_jac(
      size_t                     q,          /**< number of columns in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, columnwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, columnwise */
      )
   {
      return univariate_for_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_jac;

   /** computes sparsity of jacobian during a reverse sweep
    *
    *  For a q x 1 matrix R, we have to return the sparsity pattern of the q x 1 matrix S(x) = R * f'(x).
    *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool rev_sparse_jac(
      size_t                     q,          /**< number of rows in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, rowwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, rowwise */
      )
   {
      return univariate_rev_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<Type>::rev_sparse_hes;

   /** computes sparsity of hessian during a reverse sweep
    *
    * Assume V(x) = (g(f(x)))'' R  with f(x) = sign(x)abs(x)^p for a function g:R->R and a matrix R.
    * we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
    */
   bool rev_sparse_hes(
      const CppAD::vector<bool>& vx,         /**< indicates whether argument is a variable, or empty vector */
      const CppAD::vector<bool>& s,          /**< sparsity pattern of S = g'(y) */
      CppAD::vector<bool>&       t,          /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
      size_t                     q,          /**< number of columns in S and R */
      const CppAD::vector<bool>& r,          /**< sparsity pattern of R */
      const CppAD::vector<bool>& u,          /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
      CppAD::vector<bool>&       v           /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
      )
   {
      return univariate_rev_sparse_hes(vx, s, t, q, r, u, v);
   }

};

/** Specialization of atomic_signpower template for intervals */
template<>
class atomic_signpower<SCIPInterval> : public CppAD::atomic_base<SCIPInterval>
{
public:
   atomic_signpower()
   : CppAD::atomic_base<SCIPInterval>("signpowerint"),
     exponent(0.0)
   {
      /* indicate that we want to use bool-based sparsity pattern */
      this->option(CppAD::atomic_base<SCIPInterval>::bool_sparsity_enum);
   }

private:
   /** exponent for use in next call to forward or reverse */
   SCIP_Real exponent;

   /** stores exponent corresponding to next call to forward or reverse
    *
    *  How is this supposed to be threadsafe? (we use only one global instantiation of this class)
    */
   virtual void set_id(size_t id)
   {
      exponent = SCIPexprGetSignPowerExponent((SCIP_EXPR*)(void*)id);
   }

   /** specialization of atomic_signpower::forward template for SCIPinterval
    *
    *  @todo try to compute tighter resultants
    */
   bool forward(
      size_t                             q,  /**< lowest order Taylor coefficient that we are evaluating */
      size_t                             p,  /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<bool>&         vx, /**< indicates whether argument is a variable, or empty vector */
      CppAD::vector<bool>&               vy, /**< vector to store which function values depend on variables, or empty vector */
      const CppAD::vector<SCIPInterval>& tx, /**< values for taylor coefficients of x */
      CppAD::vector<SCIPInterval>&       ty  /**< vector to store taylor coefficients of y */
      )
   {
      assert(exponent > 0.0);
      assert(tx.size() >= p+1);
      assert(ty.size() >= p+1);
      assert(q <= p);

      if( vx.size() > 0 )
      {
         assert(vx.size() == 1);
         assert(vy.size() == 1);
         assert(p == 0);

         vy[0] = vx[0];
      }

      if( q == 0 /* q <= 0 && 0 <= p */ )
      {
         ty[0] = CppAD::signpow(tx[0], exponent);
      }

      if( q <= 1 && 1 <= p )
      {
         ty[1] = CppAD::pow(CppAD::abs(tx[0]), exponent - 1.0) * tx[1];
         ty[1] *= p;
      }

      if( q <= 2 && 2 <= p )
      {
         if( exponent != 2.0 )
         {
            ty[2]  = CppAD::signpow(tx[0], exponent - 2.0) * CppAD::square(tx[1]);
            ty[2] *= (exponent - 1.0) / 2.0;
            ty[2] += CppAD::pow(CppAD::abs(tx[0]), exponent - 1.0) * tx[2];
            ty[2] *= exponent;
         }
         else
         {
            // y'' = 2 (1/2 * sign(x) * x'^2 + |x|*x'') = sign(tx[0]) * tx[1]^2 + 2 * abs(tx[0]) * tx[2]
            ty[2]  = CppAD::sign(tx[0]) * CppAD::square(tx[1]);
            ty[2] += 2.0 * CppAD::abs(tx[0]) * tx[2];
         }
      }

      /* higher order derivatives not implemented */
      if( p > 2 )
         return false;

      return true;
   }

   /** specialization of atomic_signpower::reverse template for SCIPinterval
    *
    *  @todo try to compute tighter resultants
    */
   bool reverse(
      size_t                             p,  /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<SCIPInterval>& tx, /**< values for taylor coefficients of x */
      const CppAD::vector<SCIPInterval>& ty, /**< values for taylor coefficients of y */
      CppAD::vector<SCIPInterval>&       px, /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
      const CppAD::vector<SCIPInterval>& py  /**< values for partial derivatives of g(x) w.r.t. y */
      )
   { /*lint --e{715} */
      assert(exponent > 1);
      assert(px.size() >= p+1);
      assert(py.size() >= p+1);
      assert(tx.size() >= p+1);

      switch( p )
      {
      case 0:
         // px[0] = py[0] * p * pow(abs(tx[0]), p-1);
         px[0]  = py[0] * CppAD::pow(CppAD::abs(tx[0]), exponent - 1.0);
         px[0] *= exponent;
         break;

      case 1:
         if( exponent != 2.0 )
         {
            // px[0] = py[0] * p * abs(tx[0])^(p-1) + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
            px[0]  = py[1] * tx[1] * CppAD::signpow(tx[0], exponent - 2.0);
            px[0] *= exponent - 1.0;
            px[0] += py[0] * CppAD::pow(CppAD::abs(tx[0]), exponent - 1.0);
            px[0] *= exponent;
            // px[1] = py[1] * p * abs(tx[0])^(p-1)
            px[1]  = py[1] * CppAD::pow(CppAD::abs(tx[0]), exponent - 1.0);
            px[1] *= exponent;
         }
         else
         {
            // px[0] = py[0] * 2.0 * abs(tx[0]) + py[1] * 2.0 * sign(tx[0]) * tx[1]
            px[0]  = py[1] * tx[1] * CppAD::sign(tx[0]);
            px[0] += py[0] * CppAD::abs(tx[0]);
            px[0] *= 2.0;
            // px[1] = py[1] * 2.0 * abs(tx[0])
            px[1]  = py[1] * CppAD::abs(tx[0]);
            px[1] *= 2.0;
         }
         break;

      default:
         return false;
      }

      return true;
   }

   using CppAD::atomic_base<SCIPInterval>::for_sparse_jac;

   /** computes sparsity of jacobian during a forward sweep
    *
    * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
    * Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool for_sparse_jac(
      size_t                     q,          /**< number of columns in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, columnwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, columnwise */
      )
   {
      return univariate_for_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<SCIPInterval>::rev_sparse_jac;

   /** computes sparsity of jacobian during a reverse sweep
    *
    *  For a q x 1 matrix R, we have to return the sparsity pattern of the q x 1 matrix S(x) = R * f'(x).
    *  Since f'(x) is dense, the sparsity of S will be the sparsity of R.
    */
   bool rev_sparse_jac(
      size_t                     q,          /**< number of rows in R */
      const CppAD::vector<bool>& r,          /**< sparsity of R, rowwise */
      CppAD::vector<bool>&       s           /**< vector to store sparsity of S, rowwise */
      )
   {
      return univariate_rev_sparse_jac(q, r, s);
   }

   using CppAD::atomic_base<SCIPInterval>::rev_sparse_hes;

   /** computes sparsity of hessian during a reverse sweep
    *
    * Assume V(x) = (g(f(x)))'' R  with f(x) = sign(x)abs(x)^p for a function g:R->R and a matrix R.
    * we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
    */
   bool rev_sparse_hes(
      const CppAD::vector<bool>& vx,         /**< indicates whether argument is a variable, or empty vector */
      const CppAD::vector<bool>& s,          /**< sparsity pattern of S = g'(y) */
      CppAD::vector<bool>&       t,          /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
      size_t                     q,          /**< number of columns in S and R */
      const CppAD::vector<bool>& r,          /**< sparsity pattern of R */
      const CppAD::vector<bool>& u,          /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
      CppAD::vector<bool>&       v           /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
      )
   {
      return univariate_rev_sparse_hes(vx, s, t, q, r, u, v);
   }
};

/** template for evaluation for signpower operator */
template<class Type>
static
void evalSignPower(
   Type&                 resultant,          /**< resultant */
   const Type&           arg,                /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{
   vector<Type> in(1, arg);
   vector<Type> out(1);

   static atomic_signpower<typename Type::value_type> sp;
   sp(in, out, (size_t)(void*)expr);

   resultant = out[0];
   return;
}

#else

/** template for evaluation for signpower operator
 *
 *  Only implemented for real numbers, thus gives error by default.
 */
template<class Type>
static
void evalSignPower(
   Type&                 resultant,          /**< resultant */
   const Type&           arg,                /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{  /*lint --e{715}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalSignPower()",
      "Error: SignPower not implemented for this value type"
      );
}

/** specialization of signpower evaluation for real numbers */
template<>
void evalSignPower(
   CppAD::AD<double>&    resultant,          /**< resultant */
   const CppAD::AD<double>& arg,             /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{
   SCIP_Real exponent;

   exponent = SCIPexprGetSignPowerExponent(expr);

   if( arg == 0.0 )
      resultant = 0.0;
   else if( arg > 0.0 )
      resultant =  pow( arg, exponent);
   else
      resultant = -pow(-arg, exponent);
}

#endif


#ifndef NO_CPPAD_USER_ATOMIC

template<class Type>
SCIP_RETCODE exprEvalUser(
   SCIP_EXPR* expr,
   Type* x,
   Type& funcval,
   Type* gradient,
   Type* hessian
   )
{
   return SCIPexprEvalUser(expr, x, &funcval, gradient, hessian); /*lint !e429*/
}

template<>
SCIP_RETCODE exprEvalUser(
   SCIP_EXPR* expr,
   SCIPInterval* x,
   SCIPInterval& funcval,
   SCIPInterval* gradient,
   SCIPInterval* hessian
   )
{
   return SCIPexprEvalIntUser(expr, SCIPInterval::infinity, x, &funcval, gradient, hessian);
}

/** Automatic differentiation of user expression as CppAD user-atomic function.
 *
 * This class implements forward and reverse operations for a function given by a user expression for use within CppAD.
 */
template<class Type>
class atomic_userexpr : public CppAD::atomic_base<Type>
{
public:
   atomic_userexpr()
   : CppAD::atomic_base<Type>("userexpr"),
     expr(NULL)
   {
      /* indicate that we want to use bool-based sparsity pattern */
      this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
   }

private:
   /** user expression */
   SCIP_EXPR* expr;

   /** stores user expression corresponding to next call to forward or reverse
    *
    * how is this supposed to be threadsafe? (we use only one global instantiation of this class)
    */
   virtual void set_id(size_t id)
   {
      expr = (SCIP_EXPR*)(void*)id;
      assert(SCIPexprGetOperator(expr) == SCIP_EXPR_USER);
   }

   /** forward sweep of userexpr
    *
    * We follow http://www.coin-or.org/CppAD/Doc/atomic_forward.xml
    *   Note, that p and q are interchanged!
    *
    * For a scalar variable t, let
    *   Y(t) = f(X(t))
    *   X(t) = x^0 + x^1 t^1 + ... + x^p t^p
    * where for x^i the i an index, while for t^i the i is an exponent.
    * Thus, x^k = 1/k! X^(k) (0),   where X^(k)(.) denotes the k-th derivative.
    *
    * Next, let y^k = 1/k! Y^(k)(0) be the k'th taylor coefficient of Y. Thus,
    *   y^0 = Y^(0)(0)     =     Y(0)   = f(X(0)) = f(x^0)
    *   y^1 = Y^(1)(0)     =     Y'(0)  = f'(X(0)) * X'(0) = f'(x^0) * x^1
    *   y^2 = 1/2 Y^(2)(0) = 1/2 Y''(0) = 1/2 X'(0) * f''(X(0)) X'(0) + 1/2 * f'(X(0)) * X''(0) = 1/2 x^1 * f''(x^0) * x^1 + f'(x^0) * x^2
    *
    * As x^k = (tx[k], tx[(p+1)+k], tx[2*(p+1)+k], ..., tx[n*(p+1)+k], we get
    *   ty[0] = y^0 = f(x^0) = f(tx[{1..n}*(p+1)])
    *   ty[1] = y^1 = f'(x^0) * tx[{1..n}*(p+1)+1] = sum(i=1..n, grad[i] * tx[i*(p+1)+1]),  where grad = f'(x^0)
    *   ty[2] = 1/2 sum(i,j=1..n, x[i*(p+1)+1] * x[j*(p+1)+q] * hessian[i,j]) + sum(i=1..n, grad[i] * x[i*(p+1)+2])
    */
   bool forward(
      size_t                      q,            /**< lowest order Taylor coefficient that we are evaluating */
      size_t                      p,            /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<bool>&  vx,           /**< indicates whether argument is a variable, or empty vector */
      CppAD::vector<bool>&        vy,           /**< vector to store which function values depend on variables, or empty vector */
      const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
      CppAD::vector<Type>&        ty            /**< vector to store taylor coefficients of y */
   )
   {
      assert(expr != NULL);
      assert(ty.size() == p+1);
      assert(q <= p);

      size_t n = tx.size() / (p+1);
      assert(n == (size_t)SCIPexprGetNChildren(expr)); /*lint !e571*/
      assert(n >= 1);

      if( vx.size() > 0 )
      {
         assert(vx.size() == n);
         assert(vy.size() == 1);
         assert(p == 0);

         /* y_0 is a variable if at least one of the x_i is a variable */
         vy[0] = false;
         for( size_t i = 0; i < n; ++i )
            if( vx[i] )
            {
               vy[0] = true;
               break;
            }
      }

      Type* x = new Type[n];
      Type* gradient = NULL;
      Type* hessian = NULL;

      if( q <= 2 && 1 <= p )
         gradient = new Type[n];
      if( q <= 2 && 2 <= p )
         hessian = new Type[n*n];

      for( size_t i = 0; i < n; ++i )
         x[i] = tx[i * (p+1) + 0];  /*lint !e835*/

      if( exprEvalUser(expr, x, ty[0], gradient, hessian) != SCIP_OKAY )
      {
         delete[] x;
         delete[] gradient;
         delete[] hessian;
         return false;
      }

      if( gradient != NULL )
      {
         ty[1] = 0.0;
         for( size_t i = 0; i < n; ++i )
            ty[1] += gradient[i] * tx[i * (p+1) + 1];
      }

      if( hessian != NULL )
      {
         assert(gradient != NULL);

         ty[2] = 0.0;
         for( size_t i = 0; i < n; ++i )
         {
            for( size_t j = 0; j < n; ++j )
               ty[2] += 0.5 * hessian[i*n+j] * tx[i * (p+1) + 1] * tx[j * (p+1) + 1];

            ty[2] += gradient[i] * tx[i * (p+1) + 2];
         }
      }

      delete[] x;
      delete[] gradient;
      delete[] hessian;

      /* higher order derivatives not implemented */
      if( p > 2 )
         return false;

      return true;
   }

   /** reverse sweep of userexpr
    *
    * We follow http://www.coin-or.org/CppAD/Doc/atomic_reverse.xml
    *   Note, that there q is our p.
    *
    * For a scalar variable t, let
    *   Y(t) = f(X(t))
    *   X(t) = x^0 + x^1 t^1 + ... + x^p t^p
    * where for x^i the i an index, while for t^i the i is an exponent.
    * Thus, x^k = 1/k! X^(k) (0),   where X^(k)(.) denotes the k-th derivative.
    *
    * Next, let y^k = 1/k! Y^(k)(0) be the k'th taylor coefficient of Y. Thus,
    *   Y(t) = y^0 + y^1 t^1 + y^2 t^2 + ...
    * y^0, y^1, ... are the taylor coefficients of f(x).
    *
    * Further, let F(x^0,..,x^p) by given as F^k(x) = y^k. Thus,
    *   F^0(x) = y^0 = Y^(0)(0)   = f(x^0)
    *   F^1(x) = y^1 = Y^(1)(0)   = f'(x^0) * x^1
    *   F^2(x) = y^2 = 1/2 Y''(0) = 1/2 x^1 f''(x^0) x^1 + f'(x^0) x^2
    *
    * Given functions G: R^(p+1) -> R and H: R^(n*(p+1)) -> R, where H(x^0, x^1, .., x^p) = G(F(x^0,..,x^p)),
    * we have to return the value of \f$\partial H / \partial x^l, l = 0..p,\f$ in px. Therefor,
    * \f$
    *  px^l = \partial H / \partial x^l
    *       = sum(k=0..p, (\partial G / \partial y^k) * (\partial y^k / \partial x^l)
    *       = sum(k=0..p, py[k] * (\partial F^k / \partial x^l)
    * \f$
    *
    * For p = 0, this means
    * \f$
    *  px^0 = py[0] * \partial F^0 / \partial x^0
    *       = py[0] * \partial f(x^0) / \partial x^0
    *       = py[0] * f'(x^0)
    * \f$
    *
    * For p = 1, this means (for l = 0):
    * \f[
    *  px^0 = py[0] * \partial F^0    / \partial x^0 + py[1] * \partial F^1 / \partial x^0
    *       = py[0] * \partial f(x^0) / \partial x^0 + py[1] * \partial (f'(x^0) * x^1) / \partial x^0
    *       = py[0] * f'(x^0)                         + py[1] * f''(x^0) * x^1
    * \f]
    * and (for l=1):
    * \[
    *  px^1 = py[0] * \partial F^0    / \partial x^1 + py[1] * \partial F^1 / \partial x^1
    *       = py[0] * \partial f(x^0) / \partial x^1 + py[1] * \partial (f'(x^0) * x^1) / \partial x^0
    *       = py[0] * 0                               + py[1] * f'(x^0)
    * \f]
    *
    * As x^k = (tx[k], tx[(p+1)+k], tx[2*(p+1)+k], ..., tx[n*(p+1)+k] and
    *   px^k = (px[k], px[(p+1)+k], px[2*(p+1)+k], ..., px[n*(p+1)+k], we get
    * for p = 0:
    *   px[i] = (px^0)_i = py[0] * grad[i]
    * for p = 1:
    *   px[i*2+0] = (px^0)_i = py[0] * grad[i] + py[1] * sum(j, hessian[j,i] * tx[j*2+1])
    *   px[i*2+1] = (px^1)_i = py[1] * grad[i]
    */
   bool reverse(
      size_t                      p,            /**< highest order Taylor coefficient that we are evaluating */
      const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
      const CppAD::vector<Type>&  ty,           /**< values for taylor coefficients of y */
      CppAD::vector<Type>&        px,           /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
      const CppAD::vector<Type>&  py            /**< values for partial derivatives of g(x) w.r.t. y */
      )
   {
      assert(expr != NULL);
      assert(px.size() == tx.size());
      assert(py.size() == p+1);

      size_t n = tx.size() / (p+1);
      assert(n == (size_t)SCIPexprGetNChildren(expr)); /*lint !e571*/
      assert(n >= 1);

      Type* x = new Type[n];
      Type funcval;
      Type* gradient = new Type[n];
      Type* hessian = NULL;

      if( p == 1 )
         hessian = new Type[n*n];

      for( size_t i = 0; i < n; ++i )
         x[i] = tx[i * (p+1) + 0]; /*lint !e835*/

      if( exprEvalUser(expr, x, funcval, gradient, hessian) != SCIP_OKAY )
      {
         delete[] x;
         delete[] gradient;
         delete[] hessian;
         return false;
      }

      switch( p )
      {
      case 0:
         // px[j] = (px^0)_j = py[0] * grad[j]
         for( size_t i = 0; i < n; ++i )
            px[i] = py[0] * gradient[i];
         break;

      case 1:
         //  px[i*2+0] = (px^0)_i = py[0] * grad[i] + py[1] * sum(j, hessian[j,i] * tx[j*2+1])
         //  px[i*2+1] = (px^1)_i = py[1] * grad[i]
         assert(hessian != NULL);
         for( size_t i = 0; i < n; ++i )
         {
            px[i*2+0] = py[0] * gradient[i]; /*lint !e835*/
            for( size_t j = 0; j < n; ++j )
               px[i*2+0] += py[1] * hessian[i+n*j] * tx[j*2+1]; /*lint !e835*/

            px[i*2+1] = py[1] * gradient[i];
         }
         break;

      default:
         return false;
      }

      return true;
   } /*lint !e715*/

   using CppAD::atomic_base<Type>::for_sparse_jac;

   /** computes sparsity of jacobian during a forward sweep
    * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
    * Since we assume f'(x) to be dense, the sparsity of S will be the sparsity of R.
    */
   bool for_sparse_jac(
      size_t                     q,  /**< number of columns in R */
      const CppAD::vector<bool>& r,  /**< sparsity of R, columnwise */
      CppAD::vector<bool>&       s   /**< vector to store sparsity of S, columnwise */
      )
   {
      assert(expr != NULL);
      assert(s.size() == q);

      size_t n = r.size() / q;
      assert(n == (size_t)SCIPexprGetNChildren(expr)); /*lint !e571*/

      // sparsity for S(x) = f'(x) * R
      for( size_t j = 0; j < q; j++ )
      {
         s[j] = false;
         for( size_t i = 0; i < n; i++ )
            s[j] |= (bool)r[i * q + j]; /*lint !e1786*/
      }

      return true;
   }

   using CppAD::atomic_base<Type>::rev_sparse_jac;

   /** computes sparsity of jacobian during a reverse sweep
    * For a q x 1 matrix S, we have to return the sparsity pattern of the q x 1 matrix R(x) = S * f'(x).
    * Since we assume f'(x) to be dense, the sparsity of R will be the sparsity of S.
    */
   bool rev_sparse_jac(
      size_t                     q,  /**< number of rows in R */
      const CppAD::vector<bool>&       rt, /**< sparsity of R, rowwise */
      CppAD::vector<bool>& st  /**< vector to store sparsity of S, rowwise */
      )
   {
      assert(expr != NULL);
      assert(rt.size() == q);

      size_t n = st.size() / q;
      assert(n == (size_t)SCIPexprGetNChildren(expr));

      // sparsity for S(x)^T = f'(x)^T * R^T
      for( size_t j = 0; j < q; j++ )
         for( size_t i = 0; i < n; i++ )
            st[i * q + j] = rt[j];

      return true;
   }

   using CppAD::atomic_base<Type>::rev_sparse_hes;

   /** computes sparsity of hessian during a reverse sweep
    * Assume V(x) = (g(f(x)))'' R  for a function g:R->R and a matrix R.
    * we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
    */
   bool rev_sparse_hes(
      const CppAD::vector<bool>&              vx, /**< indicates whether argument is a variable, or empty vector */
      const CppAD::vector<bool>&              s,  /**< sparsity pattern of S = g'(y) */
      CppAD::vector<bool>&                    t,  /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
      size_t                                  q,  /**< number of columns in S and R */
      const CppAD::vector<bool>& r,  /**< sparsity pattern of R */
      const CppAD::vector<bool>& u,  /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
      CppAD::vector<bool>&      v   /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
   )
   {
      assert(expr != NULL);
      size_t n = vx.size();
      assert((size_t)SCIPexprGetNChildren(expr) == n);
      assert(s.size() == 1);
      assert(t.size() == n);
      assert(r.size() == n * q);
      assert(u.size() == q);
      assert(v.size() == n * q);

      size_t i, j, k;

      // sparsity for T(x) = S(x) * f'(x)
      for( i = 0; i < n; ++i )
         t[i] = s[0];

      // V(x) = f'(x)^T * g''(y) * f'(x) * R  +  g'(y) * f''(x) * R
      // U(x) = g''(y) * f'(x) * R
      // S(x) = g'(y)

      // back propagate the sparsity for U
      for( j = 0; j < q; j++ )
         for( i = 0; i < n; i++ )
            v[ i * q + j] = u[j];

      // include forward Jacobian sparsity in Hessian sparsity
      // sparsity for g'(y) * f''(x) * R  (Note f''(x) is assumed to be dense)
      if( s[0] )
         for( j = 0; j < q; j++ )
            for( i = 0; i < n; i++ )
               for( k = 0; k < n; ++k )
                  v[ i * q + j] |= (bool) r[ k * q + j];

      return true;
   }

};

template<class Type>
static
void evalUser(
   Type&                 resultant,          /**< resultant */
   const Type*           args,               /**< operands */
   SCIP_EXPR*            expr                /**< expression that holds the user expression */
   )
{
   assert( args != 0 );
   vector<Type> in(args, args + SCIPexprGetNChildren(expr));
   vector<Type> out(1);

   static atomic_userexpr<typename Type::value_type> u;
   u(in, out, (size_t)(void*)expr);

   resultant = out[0];
   return;
}

#else

template<class Type>
static
void evalUser(
   Type&                 resultant,          /**< resultant */
   const Type*           args,               /**< operands */
   SCIP_EXPR*            expr                /**< expression that holds the user expression */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalUser()",
      "Error: user expressions in CppAD not possible without CppAD user atomic facility"
      );
}

#endif

/** template for evaluation for minimum operator
 *
 *  Only implemented for real numbers, thus gives error by default.
 *  @todo implement own userad function
 */
template<class Type>
static
void evalMin(
   Type&                 resultant,          /**< resultant */
   const Type&           arg1,               /**< first operand */
   const Type&           arg2                /**< second operand */
   )
{  /*lint --e{715,1764}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalMin()",
      "Error: Min not implemented for this value type"
      );
}

/** specialization of minimum evaluation for real numbers */
template<>
void evalMin(
   CppAD::AD<double>&    resultant,          /**< resultant */
   const CppAD::AD<double>& arg1,            /**< first operand */
   const CppAD::AD<double>& arg2             /**< second operand */
   )
{
   resultant = MIN(arg1, arg2);
}

/** template for evaluation for maximum operator
 *
 *  Only implemented for real numbers, thus gives error by default.
 *  @todo implement own userad function
 */
template<class Type>
static
void evalMax(
   Type&                 resultant,          /**< resultant */
   const Type&           arg1,               /**< first operand */
   const Type&           arg2                /**< second operand */
   )
{  /*lint --e{715,1764}*/
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalMax()",
      "Error: Max not implemented for this value type"
      );
}

/** specialization of maximum evaluation for real numbers */
template<>
void evalMax(
   CppAD::AD<double>&    resultant,          /**< resultant */
   const CppAD::AD<double>& arg1,            /**< first operand */
   const CppAD::AD<double>& arg2             /**< second operand */
   )
{
   resultant = MAX(arg1, arg2);
}

/** template for evaluation for square-root operator
 *
 *  Default is to use the standard sqrt-function.
 */
template<class Type>
static
void evalSqrt(
   Type&                 resultant,          /**< resultant */
   const Type&           arg                 /**< operand */
   )
{
   resultant = sqrt(arg);
}

/** template for evaluation for absolute value operator */
template<class Type>
static
void evalAbs(
   Type&                 resultant,          /**< resultant */
   const Type&           arg                 /**< operand */
   )
{
   resultant = abs(arg);
}

/** specialization of absolute value evaluation for intervals
 *
 *  Use sqrt(x^2) for now @todo implement own userad function.
 */
template<>
void evalAbs(
   CppAD::AD<SCIPInterval>& resultant,       /**< resultant */
   const CppAD::AD<SCIPInterval>& arg        /**< operand */
   )
{
   vector<CppAD::AD<SCIPInterval> > in(1, arg);
   vector<CppAD::AD<SCIPInterval> > out(1);

   posintpower(in, out, 2);

   resultant = sqrt(out[0]);
}

/** integer power operation for arbitrary integer exponents */
template<class Type>
static
void evalIntPower(
   Type&                 resultant,          /**< resultant */
   const Type&           arg,                /**< operand */
   const int             exponent            /**< exponent */
   )
{
   if( exponent > 1 )
   {
      vector<Type> in(1, arg);
      vector<Type> out(1);

      posintpower(in, out, exponent);

      resultant = out[0];
      return;
   }

   if( exponent < -1 )
   {
      vector<Type> in(1, arg);
      vector<Type> out(1);

      posintpower(in, out, -exponent);

      resultant = Type(1.0)/out[0];
      return;
   }

   if( exponent == 1 )
   {
      resultant = arg;
      return;
   }

   if( exponent == 0 )
   {
      resultant = Type(1.0);
      return;
   }

   assert(exponent == -1);
   resultant = Type(1.0)/arg;
}

/** CppAD compatible evaluation of an expression for given arguments and parameters */
template<class Type>
static
SCIP_RETCODE eval(
   SCIP_EXPR*            expr,               /**< expression */
   const vector<Type>&   x,                  /**< values of variables */
   SCIP_Real*            param,              /**< values of parameters */
   Type&                 val                 /**< buffer to store expression value */
   )
{
   Type* buf = 0;

   assert(expr != NULL);

   /* todo use SCIP_MAXCHILD_ESTIMATE as in expression.c */

   if( SCIPexprGetNChildren(expr) )
   {
      if( BMSallocMemoryArray(&buf, SCIPexprGetNChildren(expr)) == NULL )  /*lint !e666*/
         return SCIP_NOMEMORY;

      for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
      {
         SCIP_CALL( eval(SCIPexprGetChildren(expr)[i], x, param, buf[i]) );
      }
   }

   switch(SCIPexprGetOperator(expr))
   {
   case SCIP_EXPR_VARIDX:
      assert(SCIPexprGetOpIndex(expr) < (int)x.size());
      val = x[SCIPexprGetOpIndex(expr)];
      break;

   case SCIP_EXPR_CONST:
      val = SCIPexprGetOpReal(expr);
      break;

   case SCIP_EXPR_PARAM:
      assert(param != NULL);
      val = param[SCIPexprGetOpIndex(expr)];
      break;

   case SCIP_EXPR_PLUS:
      assert( buf != 0 );
      val = buf[0] + buf[1];
      break;

   case SCIP_EXPR_MINUS:
      assert( buf != 0 );
      val = buf[0] - buf[1];
      break;

   case SCIP_EXPR_MUL:
      assert( buf != 0 );
      val = buf[0] * buf[1];
      break;

   case SCIP_EXPR_DIV:
      assert( buf != 0 );
      val = buf[0] / buf[1];
      break;

   case SCIP_EXPR_SQUARE:
      assert( buf != 0 );
      evalIntPower(val, buf[0], 2);
      break;

   case SCIP_EXPR_SQRT:
      assert( buf != 0 );
      evalSqrt(val, buf[0]);
      break;

   case SCIP_EXPR_REALPOWER:
      assert( buf != 0 );
      val = CppAD::pow(buf[0], SCIPexprGetRealPowerExponent(expr));
      break;

   case SCIP_EXPR_INTPOWER:
      assert( buf != 0 );
      evalIntPower(val, buf[0], SCIPexprGetIntPowerExponent(expr));
      break;

   case SCIP_EXPR_SIGNPOWER:
      assert( buf != 0 );
      evalSignPower(val, buf[0], expr);
      break;

   case SCIP_EXPR_EXP:
      assert( buf != 0 );
      val = exp(buf[0]);
      break;

   case SCIP_EXPR_LOG:
      assert( buf != 0 );
      val = log(buf[0]);
      break;

   case SCIP_EXPR_SIN:
      assert( buf != 0 );
      val = sin(buf[0]);
      break;

   case SCIP_EXPR_COS:
      assert( buf != 0 );
      val = cos(buf[0]);
      break;

   case SCIP_EXPR_TAN:
      assert( buf != 0 );
      val = tan(buf[0]);
      break;
#ifdef SCIP_DISABLED_CODE /* these operators are currently disabled */
   case SCIP_EXPR_ERF:
      assert( buf != 0 );
      val = erf(buf[0]);
      break;

   case SCIP_EXPR_ERFI:
      return SCIP_ERROR;
#endif
   case SCIP_EXPR_MIN:
      assert( buf != 0 );
      evalMin(val, buf[0], buf[1]);
      break;

   case SCIP_EXPR_MAX:
      assert( buf != 0 );
      evalMax(val, buf[0], buf[1]);
      break;

   case SCIP_EXPR_ABS:
      assert( buf != 0 );
      evalAbs(val, buf[0]);
      break;

   case SCIP_EXPR_SIGN:
      assert( buf != 0 );
      val = sign(buf[0]);
      break;

   case SCIP_EXPR_SUM:
      assert( buf != 0 );
      val = 0.0;
      for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
         val += buf[i];
      break;

   case SCIP_EXPR_PRODUCT:
      assert( buf != 0 );
      val = 1.0;
      for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
         val *= buf[i];
      break;

   case SCIP_EXPR_LINEAR:
   {
      SCIP_Real* coefs;

      coefs = SCIPexprGetLinearCoefs(expr);
      assert(coefs != NULL || SCIPexprGetNChildren(expr) == 0);

      assert( buf != 0 );
      val = SCIPexprGetLinearConstant(expr);
      for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
         val += coefs[i] * buf[i]; /*lint !e613*/
      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_Real* lincoefs;
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      SCIP_Real sqrcoef;
      Type lincoef;
      vector<Type> in(1);
      vector<Type> out(1);

      assert( buf != 0 );

      lincoefs   = SCIPexprGetQuadLinearCoefs(expr);
      nquadelems = SCIPexprGetNQuadElements(expr);
      quadelems  = SCIPexprGetQuadElements(expr);
      assert(quadelems != NULL || nquadelems == 0);

      SCIPexprSortQuadElems(expr);

      val = SCIPexprGetQuadConstant(expr);

      /* for each argument, we collect it's linear index from lincoefs, it's square coefficients and all factors from bilinear terms
       * then we compute the interval sqrcoef*x^2 + lincoef*x and add it to result */
      int i = 0;
      for( int argidx = 0; argidx < SCIPexprGetNChildren(expr); ++argidx )
      {
         if( i == nquadelems || quadelems[i].idx1 > argidx ) /*lint !e613*/
         {
            /* there are no quadratic terms with argidx in its first argument, that should be easy to handle */
            if( lincoefs != NULL )
               val += lincoefs[argidx] * buf[argidx];
            continue;
         }

         sqrcoef = 0.0;
         lincoef = lincoefs != NULL ? lincoefs[argidx] : 0.0;

         assert(i < nquadelems && quadelems[i].idx1 == argidx); /*lint !e613*/
         do
         {
            if( quadelems[i].idx2 == argidx )  /*lint !e613*/
               sqrcoef += quadelems[i].coef; /*lint !e613*/
            else
               lincoef += quadelems[i].coef * buf[quadelems[i].idx2]; /*lint !e613*/
            ++i;
         } while( i < nquadelems && quadelems[i].idx1 == argidx ); /*lint !e613*/
         assert(i == nquadelems || quadelems[i].idx1 > argidx);  /*lint !e613*/

         /* this is not as good as what we can get from SCIPintervalQuad, but easy to implement */
         if( sqrcoef != 0.0 )
         {
            in[0] = buf[argidx];
            posintpower(in, out, 2);
            val += sqrcoef * out[0];
         }

         val += lincoef * buf[argidx];
      }
      assert(i == nquadelems);

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_MONOMIAL** monomials;
      Type childval;
      Type monomialval;
      SCIP_Real exponent;
      int nmonomials;
      int nfactors;
      int* childidxs;
      SCIP_Real* exponents;
      int i;
      int j;

      assert( buf != 0 );

      val = SCIPexprGetPolynomialConstant(expr);

      nmonomials = SCIPexprGetNMonomials(expr);
      monomials  = SCIPexprGetMonomials(expr);

      for( i = 0; i < nmonomials; ++i )
      {
         nfactors  = SCIPexprGetMonomialNFactors(monomials[i]);
         childidxs = SCIPexprGetMonomialChildIndices(monomials[i]);
         exponents = SCIPexprGetMonomialExponents(monomials[i]);
         monomialval  = SCIPexprGetMonomialCoef(monomials[i]);

         for( j = 0; j < nfactors; ++j )
         {
            assert(childidxs[j] >= 0);
            assert(childidxs[j] <  SCIPexprGetNChildren(expr));

            childval = buf[childidxs[j]];
            exponent = exponents[j];

            /* cover some special exponents separately to avoid calling expensive pow function */
            if( exponent == 0.0 )
               continue;
            if( exponent == 1.0 )
            {
               monomialval *= childval;
               continue;
            }
            if( (int)exponent == exponent )
            {
               Type tmp;
               evalIntPower(tmp, childval, (int)exponent);
               monomialval *= tmp;
               continue;
            }
            if( exponent == 0.5 )
            {
               Type tmp;
               evalSqrt(tmp, childval);
               monomialval *= tmp;
               continue;
            }
            monomialval *= pow(childval, exponent);
         }

         val += monomialval;
      }

      break;
   }

   case SCIP_EXPR_USER:
      evalUser(val, buf, expr);
      break;

   case SCIP_EXPR_LAST:
   default:
      BMSfreeMemoryArrayNull(&buf);
      return SCIP_ERROR;
   }

   BMSfreeMemoryArrayNull(&buf);

   return SCIP_OKAY;
}

/** analysis an expression tree whether it requires retaping on every evaluation
 *
 *  This may be the case if the evaluation sequence depends on values of operands (e.g., in case of abs, sign, signpower, ...).
 */
static
void analyzeTree(
   SCIP_EXPRINTDATA* data,
   SCIP_EXPR*        expr
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetChildren(expr) != NULL || SCIPexprGetNChildren(expr) == 0);

   for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
      analyzeTree(data, SCIPexprGetChildren(expr)[i]);

   switch( SCIPexprGetOperator(expr) )
   {
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
#ifdef NO_CPPAD_USER_ATOMIC
   case SCIP_EXPR_SIGNPOWER:
#endif
      data->need_retape_always = true;
      break;

   case SCIP_EXPR_USER:
      data->userevalcapability &= SCIPexprGetUserEvalCapability(expr);
      break;

   default: ;
   } /*lint !e788*/

}

/** replacement for CppAD's default error handler
 *
 *  In debug mode, CppAD gives an error when an evaluation contains a nan.
 *  We do not want to stop execution in such a case, since the calling routine should check for nan's and decide what to do.
 *  Since we cannot ignore this particular error, we ignore all.
 *  @todo find a way to check whether the error corresponds to a nan and communicate this back
 */
static
void cppaderrorcallback(
   bool                  known,              /**< is the error from a known source? */
   int                   line,               /**< line where error occured */
   const char*           file,               /**< file where error occured */
   const char*           cond,               /**< error condition */
   const char*           msg                 /**< error message */
   )
{
   SCIPdebugMessage("ignore CppAD error from %sknown source %s:%d: msg: %s exp: %s\n", known ? "" : "un", file, line, msg, cond);
}

/* install our error handler */
static CppAD::ErrorHandler errorhandler(cppaderrorcallback);

/** gets name and version of expression interpreter */
const char* SCIPexprintGetName(void)
{
   return CPPAD_PACKAGE_STRING;
}

/** gets descriptive text of expression interpreter */
const char* SCIPexprintGetDesc(void)
{
   return "Algorithmic Differentiation of C++ algorithms developed by B. Bell (www.coin-or.org/CppAD)";
}

/** gets capabilities of expression interpreter (using bitflags) */
SCIP_EXPRINTCAPABILITY SCIPexprintGetCapability(
   void
   )
{
   return SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_INTFUNCVALUE |
      SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_INTGRADIENT |
      SCIP_EXPRINTCAPABILITY_HESSIAN;
}

/** creates an expression interpreter object */
SCIP_RETCODE SCIPexprintCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
   )
{
   assert(blkmem  != NULL);
   assert(exprint != NULL);

   if( BMSallocMemory(exprint) == NULL )
      return SCIP_NOMEMORY;

   (*exprint)->blkmem = blkmem;

   return SCIP_OKAY;
}

/** frees an expression interpreter object */
SCIP_RETCODE SCIPexprintFree(
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
   )
{
   assert( exprint != NULL);
   assert(*exprint != NULL);

   BMSfreeMemory(exprint);

   return SCIP_OKAY;
}

/** compiles an expression tree and stores compiled data in expression tree */
SCIP_RETCODE SCIPexprintCompile(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{ /*lint --e{429} */
   assert(tree    != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if (!data)
   {
      data = new SCIP_EXPRINTDATA();
      assert( data != NULL );
      SCIPexprtreeSetInterpreterData(tree, data);
      SCIPdebugMessage("set interpreter data in tree %p to %p\n", (void*)tree, (void*)data);
   }
   else
   {
      data->need_retape     = true;
      data->int_need_retape = true;
   }

   int n = SCIPexprtreeGetNVars(tree);

   data->X.resize(n);
   data->x.resize(n);
   data->Y.resize(1);

   data->int_X.resize(n);
   data->int_x.resize(n);
   data->int_Y.resize(1);

   if( data->root != NULL )
   {
      SCIPexprFreeDeep(exprint->blkmem, &data->root);
   }

   SCIP_EXPR* root = SCIPexprtreeGetRoot(tree);

   SCIP_CALL( SCIPexprCopyDeep(exprint->blkmem, &data->root, root) );

   data->blkmem = exprint->blkmem;

   analyzeTree(data, data->root);

   return SCIP_OKAY;
}


/** gives the capability to evaluate an expression by the expression interpreter
 *
 * In cases of user-given expressions, higher order derivatives may not be available for the user-expression,
 * even if the expression interpreter could handle these. This method allows to recognize that, e.g., the
 * Hessian for an expression is not available because it contains a user expression that does not provide
 * Hessians.
 */
SCIP_EXPRINTCAPABILITY SCIPexprintGetExprtreeCapability(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   return data->userevalcapability;
}/*lint !e715*/

/** frees interpreter data */
SCIP_RETCODE SCIPexprintFreeData(
   SCIP_EXPRINTDATA**    interpreterdata     /**< interpreter data that should freed */
   )
{
   assert( interpreterdata != NULL);
   assert(*interpreterdata != NULL);

   if( (*interpreterdata)->root != NULL )
      SCIPexprFreeDeep((*interpreterdata)->blkmem, &(*interpreterdata)->root);   

   delete *interpreterdata;
   *interpreterdata = NULL; 

   return SCIP_OKAY;
}

/** notify expression interpreter that a new parameterization is used
 *
 *  This probably causes retaping by AD algorithms.
 */
SCIP_RETCODE SCIPexprintNewParametrization(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(exprint != NULL);
   assert(tree    != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if( data != NULL )
   {
      data->need_retape     = true;
      data->int_need_retape = true;
   }

   return SCIP_OKAY;
}

/** evaluates an expression tree */
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value */
   )
{
   SCIP_EXPRINTDATA* data;

   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);
   assert(SCIPexprtreeGetNVars(tree) == (int)data->X.size());
   assert(SCIPexprtreeGetRoot(tree)  != NULL);

   int n = SCIPexprtreeGetNVars(tree);

   if( n == 0 )
   {
      SCIP_CALL( SCIPexprtreeEval(tree, NULL, val) );
      return SCIP_OKAY;
   }

   if( data->need_retape_always || data->need_retape )
   {
      for( int i = 0; i < n; ++i )
      {
         data->X[i] = varvals[i];
         data->x[i] = varvals[i];
      }

      CppAD::Independent(data->X);

      if( data->root != NULL )
         SCIP_CALL( eval(data->root, data->X, SCIPexprtreeGetParamVals(tree), data->Y[0]) );
      else
         data->Y[0] = 0.0;

      data->f.Dependent(data->X, data->Y);

      data->val = Value(data->Y[0]);
      SCIPdebugMessage("Eval retaped and computed value %g\n", data->val);

      // the following is required if the gradient shall be computed by a reverse sweep later
      // data->val = data->f.Forward(0, data->x)[0];

      data->need_retape = false;
   }
   else
   {
      assert((int)data->x.size() >= n);
      for( int i = 0; i < n; ++i )
         data->x[i] = varvals[i];

      data->val = data->f.Forward(0, data->x)[0];  /*lint !e1793*/
      SCIPdebugMessage("Eval used forward sweep to compute value %g\n", data->val);
   }

   *val = data->val;

   return SCIP_OKAY;
}

/** evaluates an expression tree on intervals */
SCIP_RETCODE SCIPexprintEvalInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables */
   SCIP_INTERVAL*        val                 /**< buffer to store interval value of expression */
   )
{
   SCIP_EXPRINTDATA* data;

   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);
   assert(SCIPexprtreeGetNVars(tree) == (int)data->int_X.size());
   assert(SCIPexprtreeGetRoot(tree)  != NULL);

   int n = SCIPexprtreeGetNVars(tree);

   if( n == 0 )
   {
      SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, NULL, val) );
      return SCIP_OKAY;
   }

   SCIPInterval::infinity = infinity;

   if( data->int_need_retape || data->need_retape_always )
   {
      for( int i = 0; i < n; ++i )
      {
         data->int_X[i] = varvals[i];
         data->int_x[i] = varvals[i];
      }

      CppAD::Independent(data->int_X);

      if( data->root != NULL )
         SCIP_CALL( eval(data->root, data->int_X, SCIPexprtreeGetParamVals(tree), data->int_Y[0]) );
      else
         data->int_Y[0] = 0.0;

      data->int_f.Dependent(data->int_X, data->int_Y);

      data->int_val = Value(data->int_Y[0]);

      data->int_need_retape = false;
   }
   else
   {
      assert((int)data->int_x.size() >= n);
      for( int i = 0; i < n; ++i )
         data->int_x[i] = varvals[i];

      data->int_val = data->int_f.Forward(0, data->int_x)[0];  /*lint !e1793*/
   }

   *val = data->int_val;

   return SCIP_OKAY;
}

/** computes value and gradient of an expression tree */
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   )
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL || new_varvals == FALSE);
   assert(val      != NULL);
   assert(gradient != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, val) );
   }
   else
      *val = data->val;

   int n = SCIPexprtreeGetNVars(tree);

   if( n == 0 )
      return SCIP_OKAY;

   vector<double> jac(data->f.Jacobian(data->x));

   for( int i = 0; i < n; ++i )
      gradient[i] = jac[i];

/* disable debug output since we have no message handler here
#ifdef SCIP_DEBUG
   SCIPdebugMessage("Grad for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("grad ="); for (int i = 0; i < n; ++i) printf("\t %g", gradient[i]); printf("\n");
#endif
*/

   return SCIP_OKAY;
}

/** computes interval value and interval gradient of an expression tree */
SCIP_RETCODE SCIPexprintGradInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable interval values changed since last call to an interval evaluation routine? */
   SCIP_INTERVAL*        val,                /**< buffer to store expression interval value */
   SCIP_INTERVAL*        gradient            /**< buffer to store expression interval gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   )
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL || new_varvals == false);
   assert(val      != NULL);
   assert(gradient != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if (new_varvals)
      SCIP_CALL( SCIPexprintEvalInt(exprint, tree, infinity, varvals, val) );
   else
      *val = data->int_val;

   int n = SCIPexprtreeGetNVars(tree);

   if( n == 0 )
      return SCIP_OKAY;

   vector<SCIPInterval> jac(data->int_f.Jacobian(data->int_x));

   for (int i = 0; i < n; ++i)
      gradient[i] = jac[i];

/* disable debug output since we have no message handler here
#ifdef SCIP_DEBUG
   SCIPdebugMessage("GradInt for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t [%g,%g]", SCIPintervalGetInf(data->int_x[i]), SCIPintervalGetSup(data->int_x[i])); printf("\n");
   SCIPdebugMessage("grad ="); for (int i = 0; i < n; ++i) printf("\t [%g,%g]", SCIPintervalGetInf(gradient[i]), SCIPintervalGetSup(gradient[i])); printf("\n");
#endif
*/

   return SCIP_OKAY;
}

/** gives sparsity pattern of hessian
 *
 *  NOTE: this function might be replaced later by something nicer.
 *  Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Bool*            sparsity            /**< buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */
   )
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL);
   assert(sparsity != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   int n = SCIPexprtreeGetNVars(tree);
   if( n == 0 )
      return SCIP_OKAY;

   int nn = n*n;

   if( data->need_retape_always )
   {
      // @todo can we do something better here, e.g., by looking at the expression tree by ourself?

      for( int i = 0; i < nn; ++i )
         sparsity[i] = TRUE;

/* disable debug output since we have no message handler here
#ifdef SCIP_DEBUG
      SCIPdebugMessage("HessianSparsityDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
      SCIPdebugMessage("sparsity = all elements, due to discontinuouities\n");
#endif
*/

      return SCIP_OKAY;
   }

   if( data->need_retape )
   {
      SCIP_Real val;
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, &val) );
   }

   SCIPdebugMessage("calling ForSparseJac\n");

   vector<bool> r(nn, false);
   for (int i = 0; i < n; ++i)
      r[i*n+i] = true;  /*lint !e647 !e1793*/
   (void) data->f.ForSparseJac(n, r); // need to compute sparsity for Jacobian first

   SCIPdebugMessage("calling RevSparseHes\n");

   vector<bool> s(1, true);
   vector<bool> sparsehes(data->f.RevSparseHes(n, s));

   for( int i = 0; i < nn; ++i )
      sparsity[i] = sparsehes[i];

/* disable debug output since we have no message handler here
#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianSparsityDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("sparsity ="); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (sparsity[i*n+j]) printf(" (%d,%d)", i, j); printf("\n");
#endif
*/

   return SCIP_OKAY;
}

/** computes value and dense hessian of an expression tree
 *
 *  The full hessian is computed (lower left and upper right triangle).
 */
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            hessian             /**< buffer to store hessian values, need to have size at least n*n */
   )
{
   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL || new_varvals == FALSE);
   assert(val     != NULL);
   assert(hessian != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, val) );
   }
   else
      *val = data->val;

   int n = SCIPexprtreeGetNVars(tree);

   if( n == 0 )
      return SCIP_OKAY;

#if 1
   /* this one uses reverse mode */
   vector<double> hess(data->f.Hessian(data->x, 0));

   int nn = n*n;
   for (int i = 0; i < nn; ++i)
      hessian[i] = hess[i];

#else
   /* this one uses forward mode */
   for( int i = 0; i < n; ++i )
      for( int j = 0; j < n; ++j )
      {
         vector<int> ii(1,i);
         vector<int> jj(1,j);
         hessian[i*n+j] = data->f.ForTwo(data->x, ii, jj)[0];
      }
#endif

/* disable debug output since we have no message handler here
#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("hess ="); for (int i = 0; i < n*n; ++i) printf("\t %g", hessian[i]); printf("\n");
#endif
*/
   return SCIP_OKAY;
}
