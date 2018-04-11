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

/**@file   lpi_spx2.cpp
 * @ingroup LPIS
 * @brief  LP interface for SoPlex version 2.0 and higher
 * @author Matthias Miltenberger
 * @author Ambros Gleixner
 *
 * This is an implementation of SCIP's LP interface for SoPlex using the extended and improved interface of SoPlex 2.0
 *
 * For debugging purposes, the SoPlex results can be double checked with CPLEX if WITH_LPSCHECK is defined. This may
 * yield false positives, since the LP is dumped to a file for transfering it to CPLEX, hence, precision may be lost.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define STRONGBRANCH_RESTOREBASIS            /**< if defined then in SCIPlpiStrongbranch() we restore the basis after the
                                              *   down branch and after the up branch; if false only after the end of a
                                              *   strong branching phase, which however seems to mostly increase strong
                                              *   branching time and iterations */

/* in this case the SoPlex results are double checked using CPLEX */
#ifdef WITH_LPSCHECK
#include <cplex.h>

#define CHECK_SPXSOLVE                  true /**< shall the SoPlex results in spxSolve() be double checked using CPLEX? */
#define CHECK_SPXSTRONGBRANCH           true /**< shall the SoPlex results in SCIPlpStrongbranch() be double checked using CPLEX? */
#define CHECK_START                     0    /**< skip first CHECK_START number of checks */
#define EXIT_AT_WRONG_RESULT            false/**< shall program be exited if CPLEX returns different result than SoPlex? */
#define EXIT_AT_CPXERROR                false/**< shall program be exited if CPLEX returns an error? */

#define CPX_CALL(x)                     do                                                                                  \
                                        {                                                                                   \
                                           int _cpxstat_;               \
                                           if( (_cpxstat_ = (x)) != 0 )                                                     \
                                           {                                                                                \
                                              SCIPmessagePrintWarning(_messagehdlr, "CPLEX error <%d>; SoPlex result unchecked\n", _cpxstat_); \
                                              if( EXIT_AT_CPXERROR )                                                        \
                                              {                                                                             \
                                                 exit(1);                                                                   \
                                              }                                                                             \
                                              else                                                                          \
                                              {                                                                             \
                                                 goto ENDCHECK;                                                             \
                                              }                                                                             \
                                           }                                                                                \
                                        }                                                                                   \
                                        while( false )
#endif

/* check the return value of setParam methods */
#define CHECK_SOPLEX_PARAM(x)                                                           \
   if( !x )                                                                             \
   {                                                                                    \
      SCIPmessagePrintWarning(_messagehdlr, "SoPlex: unsupported parameter value\n");   \
   }

/* remember the original value of the SCIP_DEBUG define and undefine it */
#ifdef SCIP_DEBUG
#define ___DEBUG
#undef SCIP_DEBUG
#endif

/* include SoPlex solver */
#include "soplex.h"

/* define subversion for versions <= 1.5.0.1 */
#ifndef SOPLEX_SUBVERSION
#define SOPLEX_SUBVERSION 0
#endif
/* define API version for versions <= 3.0.0 */
#ifndef SOPLEX_APIVERSION
#define SOPLEX_APIVERSION 0
#endif

/* check version */
#if (SOPLEX_VERSION < 200 || (SOPLEX_VERSION == 200 && SOPLEX_SUBVERSION < 2) || (SOPLEX_VERSION > 200 && SOPLEX_VERSION < 201))
#error "This interface is not compatible with SoPlex versions prior to 2.0.0.2"
#endif

#include "spxgithash.h"

/* reset the SCIP_DEBUG define to its original SCIP value */
#undef SCIP_DEBUG
#ifdef ___DEBUG
#define SCIP_DEBUG
#undef ___DEBUG
#endif

/* define snprintf when using a too old MSVC version */
#if defined(_MSC_VER) && _MSC_VER < 1900
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif

#define SOPLEX_VERBLEVEL                5    /**< verbosity level for LPINFO */

#include "scip/pub_message.h"

/********************************************************************/
/*----------------------------- C++ --------------------------------*/
/********************************************************************/

/* in C++ we have to use "0" instead of "(void*)0" */
#undef NULL
#define NULL 0

#include <cassert>
using namespace soplex;


/** Macro for a single SoPlex call for which exceptions have to be catched - return an LP error. We
 *  make no distinction between different exception types, e.g., between memory allocation and other
 *  exceptions.
 */
#ifndef NDEBUG
#define SOPLEX_TRY(messagehdlr, x)  do                                  \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch( const SPxMemoryException& E )                              \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPerrorMessage("SoPlex threw a memory exception: %s\n", s.c_str()); \
         return SCIP_ERROR;                                             \
      }                                                                 \
      catch( const SPxException& E )                                    \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPmessagePrintWarning((messagehdlr), "SoPlex threw an exception: %s\n", s.c_str()); \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }                                                                    \
   while( FALSE )

#else
#define SOPLEX_TRY(messagehdlr, x)  do                                  \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch( const SPxMemoryException& E )                              \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPerrorMessage("SoPlex threw a memory exception: %s\n", s.c_str()); \
         return SCIP_ERROR;                                             \
      }                                                                 \
      catch( const SPxException& )                                      \
      {                                                                 \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }                                                                    \
   while( FALSE )
#endif

/* Macro for a single SoPlex call for which exceptions have to be catched - abort if they
 * arise. SCIP_ABORT() is not accessible here.
 */
#define SOPLEX_TRY_ABORT(x)  do                                         \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch( const SPxException& E )                                    \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPerrorMessage("SoPlex threw an exception: %s\n", s.c_str()); \
         abort();                                                       \
      }                                                                 \
   }                                                                    \
   while( FALSE )



/** SCIP's SoPlex class */
class SPxSCIP : public SoPlex
{/*lint !e1790*/
   bool                  _lpinfo;
   bool                  _fromscratch;
   char*                 _probname;
   DataArray<SPxSolver::VarStatus> _colStat;  /**< column basis status used for strong branching */
   DataArray<SPxSolver::VarStatus> _rowStat;  /**< row basis status used for strong branching */
#ifdef WITH_LPSCHECK
   int                   _checknum;
   bool                  _doublecheck;
   CPXENVptr             _cpxenv;           /**< CPLEX memory environment */
   CPXLPptr              _cpxlp;            /**< CPLEX lp structure */
#endif
   SCIP_MESSAGEHDLR*     _messagehdlr;      /**< messagehdlr handler for printing messages, or NULL */

public:
   SPxSCIP(
      SCIP_MESSAGEHDLR*  messagehdlr = NULL, /**< message handler */
      const char*        probname = NULL     /**< name of problem */
      )
      : _lpinfo(false),
        _fromscratch(false),
        _probname(NULL),
        _colStat(0),
        _rowStat(0),
        _messagehdlr(messagehdlr)
   {
      if ( probname != NULL )
         SOPLEX_TRY_ABORT( setProbname(probname) );

#ifdef WITH_LPSCHECK
      int cpxstat;
      _checknum = 0;
      _doublecheck = false;
      _cpxenv = CPXopenCPLEX(&cpxstat);
      assert(_cpxenv != NULL);
      _cpxlp = CPXcreateprob(_cpxenv, &cpxstat, probname != NULL ? probname : "spxcheck");
      (void) CPXsetintparam(_cpxenv, CPX_PARAM_SCRIND, 0);
#endif
   }

   virtual ~SPxSCIP()
   {
      if( _probname != NULL )
         spx_free(_probname); /*lint !e1551*/

      freePreStrongbranchingBasis(); /*lint !e1551*/

#ifdef WITH_LPSCHECK
      (void) CPXfreeprob(_cpxenv, &_cpxlp);
      (void) CPXcloseCPLEX(&_cpxenv);
#endif
   }/*lint -e1579*/

   // we might need these methods to return the original values SCIP provided, even if they could not be set
   /** return feastol set by SCIPlpiSetRealpar(), which might be tighter than what SoPlex accepted */
   Real feastol() const
   {
      return realParam(FEASTOL);
   }

   /** set feastol and store value in case SoPlex only accepts a larger tolerance */
   void setFeastol(
      const Real d
      )
   {
      CHECK_SOPLEX_PARAM(setRealParam(FEASTOL, d));
   }

   /** return opttol set by SCIPlpiSetRealpar(), which might be tighter than what SoPlex accepted */
   Real opttol() const
   {
      return realParam(OPTTOL);
   }

   /** set opttol and store value in case SoPlex only accepts a larger tolerance */
   void setOpttol(
      const Real d
      )
   {
      CHECK_SOPLEX_PARAM(setRealParam(OPTTOL, d));
   }

   /** get objective limit according to objective sense */
   Real getObjLimit() const
   {
      return (intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
         ? realParam(SoPlex::OBJLIMIT_UPPER)
         : realParam(SoPlex::OBJLIMIT_LOWER);
   }

   // @todo realize this with a member variable as before
   bool getFromScratch() const
   {
      return _fromscratch;
   }

   void setFromScratch(bool fs)
   {
      _fromscratch = fs;
   }

   // @todo member variable?
   bool getLpInfo() const
   {
      return _lpinfo;
   }

   void setLpInfo(bool lpinfo)
   {
      _lpinfo = lpinfo;
   }

   // @todo member variable?
   void setProbname(const char* probname)
   {
      int len;

      assert(probname != NULL);
      if( _probname != NULL )
         spx_free(_probname);
      len = (int)strlen(probname);
      spx_alloc(_probname, len + 1);
      strncpy(_probname, probname, len); /*lint !e732*/
      _probname[len] = '\0';
   }

   void setRep(SPxSolver::Representation p_rep)
   {
      if( p_rep == SPxSolver::COLUMN && intParam(REPRESENTATION) == REPRESENTATION_ROW )
      {
         SCIPdebugMessage("switching to column representation of the basis\n");
         CHECK_SOPLEX_PARAM(setIntParam(REPRESENTATION, REPRESENTATION_COLUMN));
      }
      else if( (p_rep == SPxSolver::ROW && intParam(REPRESENTATION) == REPRESENTATION_COLUMN) )
      {
         SCIPdebugMessage("switching to row representation of the basis\n");
         CHECK_SOPLEX_PARAM(setIntParam(REPRESENTATION, REPRESENTATION_ROW));
      }
   }

#ifdef WITH_LPSCHECK
   bool getDoubleCheck()
   {
      _checknum++;
      return _doublecheck && _checknum + 1 >= CHECK_START;
   }

   void setDoubleCheck(bool dc)
   {
      _doublecheck = dc;
   }

   const char* spxStatusString(const SPxSolver::Status stat) const
   {
      switch( stat )
      {
      case SPxSolver::ABORT_TIME:
         return "ABORT_TIME";
      case SPxSolver::ABORT_ITER:
         return "ABORT_ITER";
      case SPxSolver::ABORT_VALUE:
         return "ABORT_VALUE";
      case SPxSolver::SINGULAR:
         return "SINGULAR";
      case SPxSolver::REGULAR:
         return "REGULAR";
      case SPxSolver::UNKNOWN:
         return "UNKNOWN";
      case SPxSolver::OPTIMAL:
         return "OPTIMAL";
      case SPxSolver::UNBOUNDED:
         return "UNBOUNDED";
      case SPxSolver::INFEASIBLE:
         return "INFEASIBLE";
      default:
         return "UNKNOWN";
      }  /*lint !e788*/
   }

   const char* cpxStatusString(const int stat) const
   {
      switch( stat )
      {
      case CPX_STAT_ABORT_TIME_LIM:
         return "ABORT_TIME";
      case CPX_STAT_ABORT_IT_LIM:
         return "ABORT_ITER";
      case CPX_STAT_ABORT_OBJ_LIM:
         return "ABORT_VALUE";
      case CPX_STAT_OPTIMAL:
         return "OPTIMAL";
      case CPX_STAT_OPTIMAL_INFEAS:
         return "CPX_STAT_OPTIMAL_INFEAS: OPT SOL INFEASIBLE AFTER UNSCALING";
      case CPX_STAT_UNBOUNDED:
         return "UNBOUNDED";
      case CPX_STAT_INFEASIBLE:
         return "INFEASIBLE";
      case CPX_STAT_INForUNBD:
         return "INFEASIBLE or UNBOUNDED";
      case CPX_STAT_NUM_BEST:
         return "CPX_STAT_NUM_BEST: SOL AVAILABLE BUT NOT PROVEN OPTIMAL DUE TO NUM TROUBLE";
      default:
         return "UNKNOWN";
      }  /*lint !e788*/
   }
#endif

#ifndef NDEBUG
   bool checkConsistentBounds() const
   {
      for( int i = 0; i < numColsReal(); ++i )
      {
         if( lowerReal(i) > upperReal(i) + realParam(SoPlex::EPSILON_ZERO) )
         {
            SCIPerrorMessage("inconsistent bounds on column %d: lower=%.17g, upper=%.17g\n",
               i, lowerReal(i), upperReal(i));
            return false;
         }
      }

      return true;
   }

   bool checkConsistentSides() const
   {
      for( int i = 0; i < numRowsReal(); ++i )
      {
         if( lhsReal(i) > rhsReal(i) + realParam(SoPlex::EPSILON_ZERO) )
         {
            SCIPerrorMessage("inconsistent sides on row %d: lhs=%.17g, rhs=%.17g\n",
               i, lhsReal(i), rhsReal(i));
            return false;
         }
      }

      return true;
   }
#endif

   void trySolve(bool printwarning = true)
   {
      Real timespent;
      Real timelimit;

      try
      {
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
         (void) optimize();
#else
         (void) solve();
#endif
      }
      catch(const SPxException& x)
      {
         std::string s = x.what();
         if( printwarning )
         {
            SCIPmessagePrintWarning(_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
         }

         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(status() != SPxSolver::OPTIMAL);
      }

      assert(intParam(ITERLIMIT) < 0 || numIterations() <= intParam(ITERLIMIT));

      /* update time limit */
      timespent = solveTime();
      if( timespent > 0 )
      {
         /* get current time limit */
         timelimit = realParam(TIMELIMIT);
         if( timelimit > timespent )
            timelimit -= timespent;
         else
            timelimit = 0;
         /* set new time limit */
         assert(timelimit >= 0);
         CHECK_SOPLEX_PARAM(setRealParam(TIMELIMIT, timelimit));
      }
   }

   SPxSolver::Status doSolve(bool printwarning = true)
   {
      SPxOut::Verbosity verbosity;

      SPxSolver::Status spxStatus;

      /* store and set verbosity */
      verbosity = spxout.getVerbosity();
      spxout.setVerbosity((SPxOut::Verbosity)(getLpInfo() ? SOPLEX_VERBLEVEL : 0));

      assert(checkConsistentBounds());
      assert(checkConsistentSides());

#ifdef WITH_LPSCHECK
      /* dump LP with current basis and settings saved in SoPlex */
      if( getDoubleCheck() )
         writeStateReal("spxcheck", NULL, NULL);
#endif

      trySolve(printwarning);
      spxStatus = status();

      /* for safety reset iteration limit */
//       setTerminationIter(_itlim);

#ifdef WITH_LPSCHECK
      bool minimize = intParam(OBJSENSE) == OBJSENSE_MINIMIZE;
      Real objLimitUpper = realParam(OBJLIMIT_UPPER);
      Real objLimitLower = realParam(OBJLIMIT_LOWER);

      /* if SoPlex gave a definite answer, we double check if it is consistent with CPLEX's answer */
      if( getDoubleCheck() && (spxStatus == SPxSolver::OPTIMAL || spxStatus == SPxSolver::UNBOUNDED || spxStatus == SPxSolver::INFEASIBLE || spxStatus == SPxSolver::ABORT_VALUE) )
      {
         SCIP_Real cpxobj;
         int cpxstat;

         /* read LP with basis */
         CPX_CALL( CPXreadcopyprob(_cpxenv, _cpxlp, "spxcheck.mps", NULL) );
         CPX_CALL( CPXreadcopybase(_cpxenv, _cpxlp, "spxcheck.bas") );

         /* set tolerances */
         CPX_CALL( CPXsetdblparam(_cpxenv, CPX_PARAM_EPOPT, MAX(opttol(), 1e-9)) );
         CPX_CALL( CPXsetdblparam(_cpxenv, CPX_PARAM_EPRHS, MAX(feastol(), 1e-9)) );

         /* solve LP */
         CPX_CALL( CPXlpopt(_cpxenv, _cpxlp) );

         /* get solution status and objective value */
         CPX_CALL( CPXsolution(_cpxenv, _cpxlp, &cpxstat, &cpxobj, NULL, NULL, NULL, NULL) );
         if( !minimize )
            cpxobj *= -1.0;

         /* check for inconsistent statuses */
         if( cpxstat == CPX_STAT_OPTIMAL_INFEAS )
         {
            SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s)\n",
               _probname, spxStatus, spxStatusString(spxStatus), cpxstat, cpxStatusString(cpxstat));
            if( EXIT_AT_CPXERROR )
               exit(1);
         }
         else if( (spxStatus == SPxSolver::OPTIMAL && cpxstat != CPX_STAT_OPTIMAL)
            || (spxStatus == SPxSolver::UNBOUNDED && cpxstat != CPX_STAT_UNBOUNDED)
            || (spxStatus == SPxSolver::INFEASIBLE && cpxstat != CPX_STAT_INFEASIBLE) )
         {
            SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s) (checknum=%d)\n",
               _probname, spxStatus, spxStatusString(spxStatus), cpxstat, cpxStatusString(cpxstat), _checknum);
            if( EXIT_AT_WRONG_RESULT )
               exit(1);
         }
         else if( spxStatus == SPxSolver::ABORT_VALUE )
         {
            switch( cpxstat )
            {
            case CPX_STAT_OPTIMAL:
               if( (minimize && LTrel(cpxobj, objLimitUpper, 2*opttol()))
                  || (!minimize && GTrel(cpxobj, objLimitLower, 2*opttol())) )
               {
                  SCIPerrorMessage("In %s: SoPlex returned status=%d (%s) while CPLEX claims obj=%.10f %s %.10f=obj.limit (%s) (checknum=%d)\n",
                     _probname, spxStatus, spxStatusString(spxStatus), cpxobj, minimize ? "<" : ">",
                     minimize ? objLimitUpper : objLimitLower, cpxStatusString(cpxstat), _checknum);
                  if( EXIT_AT_WRONG_RESULT )
                     exit(1);
               }
               else if( (minimize && cpxobj < objLimitUpper) || (!minimize && cpxobj > objLimitLower) )
               {
                  SCIPerrorMessage("In %s: SoPlex returned status=%d (%s) while CPLEX claims obj=%.10f %s %.10f=obj.limit (%s) (checknum=%d)\n",
                     _probname, spxStatus, spxStatusString(spxStatus), cpxobj, minimize? "<" : ">",
                     minimize ? objLimitUpper : objLimitLower, cpxStatusString(cpxstat), _checknum);
               }
               break;
            case CPX_STAT_OPTIMAL_INFEAS:
            case CPX_STAT_NUM_BEST:
               if( (minimize && cpxobj < objLimitUpper) || (!minimize && cpxobj > objLimitLower) )
               {
                  SCIPerrorMessage("In %s: SoPlex returned status=%d (%s) while CPLEX claims obj=%.10f %s %.10f=obj.limit (%s) (checknum=%d)\n",
                     _probname, spxStatus, spxStatusString(spxStatus), cpxobj, minimize ? "<" : ">",
                     minimize ? objLimitUpper : objLimitLower, cpxStatusString(cpxstat), _checknum);
               }
               break;
            case CPX_STAT_INFEASIBLE:
               break;
            case CPX_STAT_UNBOUNDED:
               SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s) (checknum=%d)\n",
                  _probname, spxStatus, spxStatusString(spxStatus), cpxstat, cpxStatusString(cpxstat), _checknum);
               if( EXIT_AT_WRONG_RESULT )
                  exit(1);
               break;
            case CPX_STAT_INForUNBD:
            default:
               SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s) (checknum=%d)\n",
                  _probname, spxStatus, spxStatusString(spxStatus), cpxstat, cpxStatusString(cpxstat), _checknum);
               break;
            }  /*lint !e788*/
         }
         /* check for same objective values */
         else if( spxStatus == SPxSolver::OPTIMAL )
         {
            if( (minimize && LTrel(objValueReal(), cpxobj, 2*opttol()))
               || (!minimize && GTrel(objValueReal(), cpxobj, 2*opttol())) )
            {
               /* SCIPerrorMessage("In %s: LP optimal; SoPlex value=%.10f %s CPLEX value=%.10f too good (checknum=%d)\n", value(),
                  _probname, getSense() == SPxSolver::MINIMIZE ? "<" : ">", cpxobj, _checknum); */
            }
            else if( (minimize && GTrel(objValueReal(), cpxobj, 2*opttol()))
               || (!minimize && LTrel(objValueReal(), cpxobj, 2*opttol())) )
            {
               SCIPerrorMessage("In %s: LP optimal; SoPlex value=%.10f %s CPLEX value=%.10f suboptimal (checknum=%d)\n", objValueReal(),
                  _probname, minimize ? ">" : "<", cpxobj, _checknum);
               if( EXIT_AT_WRONG_RESULT )
                  exit(1);
            }
         }
      }

   ENDCHECK:
#endif

      /* restore verbosity */
      spxout.setVerbosity(verbosity);

      return spxStatus;
   }

   /** save the current basis */
   void savePreStrongbranchingBasis()
   {
      _rowStat.reSize(numRowsReal());
      _colStat.reSize(numColsReal());

      try
      {
         getBasis(_rowStat.get_ptr(), _colStat.get_ptr());
      }
#ifndef NDEBUG
      catch(const SPxException& x)
      {
         std::string s = x.what();
         SCIPmessagePrintWarning(_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());

         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(status() != SPxSolver::OPTIMAL);
      }
#else
      catch(const SPxException&)
      { }
#endif
   }

   /** restore basis */
   void restorePreStrongbranchingBasis()
   {
      assert(_rowStat.size() == numRowsReal());
      assert(_colStat.size() == numColsReal());

      try
      {
         setBasis(_rowStat.get_ptr(), _colStat.get_ptr());
      }
#ifndef NDEBUG
      catch(const SPxException& x)
      {
         std::string s = x.what();
         SCIPmessagePrintWarning(_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
      catch(const SPxException&)
      {
#endif
         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(status() != SPxSolver::OPTIMAL);
      }
   }

   /** if basis is in store, delete it without restoring it */
   void freePreStrongbranchingBasis()
   {
      _rowStat.clear();
      _colStat.clear();
   }

   /** is pre-strong-branching basis freed? */
   bool preStrongbranchingBasisFreed() const
   {
      return ((_rowStat.size() == 0 ) && (_colStat.size() == 0));
   }

   /** provides access for temporary storage of basis status of rows */
   DataArray<SPxSolver::VarStatus>& rowStat()
   {
      return _rowStat; /*lint !e1536*/
   }

   /** provides access for temporary storage of basis status or columns */
   DataArray<SPxSolver::VarStatus>& colStat()
   {
      return _colStat; /*lint !e1536*/
   }

}; /*lint !e1748*/




/********************************************************************/
/*-----------------------------  C  --------------------------------*/
/********************************************************************/

#include "lpi/lpi.h"
#include "scip/bitencode.h"

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE



/** LP interface */
struct SCIP_LPi
{
   SPxSCIP*              spx;                /**< our SoPlex implementation */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   SCIP_PRICING          pricing;            /**< current pricing strategy */
   SCIP_Bool             solved;             /**< was the current LP solved? */
   SCIP_Real             conditionlimit;     /**< maximum condition number of LP basis counted as stable (-1.0: no limit) */
   SCIP_Bool             checkcondition;     /**< should condition number of LP basis be checked for stability? */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};

/** LPi norms to store dual steepest edge */
struct SCIP_LPiNorms
{
   int                   nrows;              /**< number of stored norms corresponding to rows */
   int                   ncols;              /**< number of stored norms corresponding to cols */
   SCIP_Real*            norms;              /**< norms to be (re)stored */
};



/*
 * dynamic memory arrays
 */

/** resizes cstat array to have at least num entries */
static
SCIP_RETCODE ensureCstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->cstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->cstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->cstat, newsize) );
      lpi->cstatsize = newsize;
   }
   assert(num <= lpi->cstatsize);

   return SCIP_OKAY;
}

/** resizes rstat array to have at least num entries */
static
SCIP_RETCODE ensureRstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->rstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->rstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rstat, newsize) );
      lpi->rstatsize = newsize;
   }
   assert(num <= lpi->rstatsize);

   return SCIP_OKAY;
}




/*
 * LPi state methods
 */

/** returns the number of packets needed to store column packet information */
static
int colpacketNum(
   int                   ncols               /**< number of columns to store */
   )
{
   return (ncols+(int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows+(int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET;
}

/** store row and column basis status in a packed LPi state object */
static
void lpistatePack(
   SCIP_LPISTATE*       lpistate,            /**< pointer to LPi state data */
   const int*           cstat,               /**< basis status of columns in unpacked format */
   const int*           rstat                /**< basis status of rows in unpacked format */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPencodeDualBit(cstat, lpistate->packcstat, lpistate->ncols);
   SCIPencodeDualBit(rstat, lpistate->packrstat, lpistate->nrows);
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   const SCIP_LPISTATE* lpistate,            /**< pointer to LPi state data */
   int*                 cstat,               /**< buffer for storing basis status of columns in unpacked format */
   int*                 rstat                /**< buffer for storing basis status of rows in unpacked format */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPdecodeDualBit(lpistate->packcstat, cstat, lpistate->ncols);
   SCIPdecodeDualBit(lpistate->packrstat, rstat, lpistate->nrows);
}

/** creates LPi state information object */
static
SCIP_RETCODE lpistateCreate(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncols,              /**< number of columns to store */
   int                   nrows               /**< number of rows to store */
   )
{
   assert(lpistate != NULL);
   assert(blkmem != NULL);
   assert(ncols >= 0);
   assert(nrows >= 0);

   int nColPackets = colpacketNum(ncols);
   int nRowPackets = rowpacketNum(nrows);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packcstat, nColPackets) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packrstat, nRowPackets) );

   return SCIP_OKAY;
}

/** frees LPi state information */
static
void lpistateFree(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state information (like basis information) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);
   assert(lpistate != NULL);
   assert(*lpistate != NULL);

   int nColPackets = colpacketNum((*lpistate)->ncols);
   int nRowPackets = rowpacketNum((*lpistate)->nrows);

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, nColPackets);
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, nRowPackets);
   BMSfreeBlockMemory(blkmem, lpistate);
}




/*
 * local methods
 */


/** marks the current LP to be unsolved */
static
void invalidateSolution(SCIP_LPI* lpi)
{
   assert(lpi != NULL);
   lpi->solved = FALSE;
}



/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static char spxname[100];
static char spxdesc[200];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolverName()\n");

#if (SOPLEX_SUBVERSION > 0)
   snprintf(spxname, 100, "SoPlex %d.%d.%d.%d", SOPLEX_VERSION/100, (SOPLEX_VERSION % 100)/10, SOPLEX_VERSION % 10, SOPLEX_SUBVERSION); /*lint !e778 !e845*/
#else
   snprintf(spxname, 100, "SoPlex %d.%d.%d", SOPLEX_VERSION/100, (SOPLEX_VERSION % 100)/10, SOPLEX_VERSION % 10); /*lint !e778 !e845*/
#endif
   return spxname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   snprintf(spxdesc, 200, "%s [GitHash: %s]", "Linear Programming Solver developed at Zuse Institute Berlin (soplex.zib.de)"
#ifdef WITH_LPSCHECK
     " - including CPLEX double check"
#endif
   , getGitHash());

   return spxdesc;
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->spx;
}

/** pass integrality information about variables to the solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0.  */
   )
{
#if (SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 3))
   assert(ncols == lpi->spx->numColsReal() || (ncols == 0 && intInfo == NULL));
   lpi->spx->setIntegralityInformation(ncols, intInfo);
   return SCIP_OKAY;
#else
   SCIPerrorMessage("SCIPlpiSetIntegralityInformation() has not been implemented yet.\n");
   return SCIP_LPERROR;
#endif
}

/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert(lpi != NULL);
   assert(name != NULL);

   /* create SoPlex object */
   SCIP_ALLOC( BMSallocMemory(lpi) );

   /* we use this construction to allocate the memory for the SoPlex class also via the blockmemshell */
   (*lpi)->spx = static_cast<SPxSCIP*>(BMSallocMemoryCPP(sizeof(SPxSCIP)));
   SOPLEX_TRY( messagehdlr, (*lpi)->spx = new ((*lpi)->spx) SPxSCIP(messagehdlr, name) );
   (void) (*lpi)->spx->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_ONLYREAL);
   (void) (*lpi)->spx->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_REAL);
   (void) (*lpi)->spx->setIntParam(SoPlex::REPRESENTATION, SoPlex::REPRESENTATION_AUTO);

   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->conditionlimit = -1.0;
   (*lpi)->checkcondition = FALSE;
   (*lpi)->messagehdlr = messagehdlr;

   invalidateSolution(*lpi);

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   /* set default pricing */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRICING, (int)(*lpi)->pricing) );

   {
      SPxOut::Verbosity verbosity = (*lpi)->spx->spxout.getVerbosity();
      (*lpi)->spx->spxout.setVerbosity((SPxOut::Verbosity)((*lpi)->spx->getLpInfo() ? SOPLEX_VERBLEVEL : 0));
      (*lpi)->spx->printVersion();
      (*lpi)->spx->spxout.setVerbosity(verbosity);
   }

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert((*lpi)->spx != NULL);

   /* free LP using destructor and free memory via blockmemshell */
   (*lpi)->spx->~SPxSCIP();
   BMSfreeMemory(&((*lpi)->spx));

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                /*rownames*/,       /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
#ifndef NDEBUG
   {
      int j;
      for( j = 0; j < nnonz; j++ )
         assert( val[j] != 0 );
   }
#endif

   SCIPdebugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(beg != NULL);
   assert(ind != NULL);
   assert(val != NULL);

   invalidateSolution(lpi);
   assert(lpi->spx->preStrongbranchingBasisFreed());

   try
   {
      SPxSCIP* spx = lpi->spx;
      LPRowSet rows(nrows);
      DSVector emptyVector(0);
      int i;

      spx->clearLPReal();

      /* set objective sense */
      (void) spx->setIntParam(SoPlex::OBJSENSE, (objsen == SCIP_OBJSEN_MINIMIZE ? SoPlex::OBJSENSE_MINIMIZE : SoPlex::OBJSENSE_MAXIMIZE));

      /* create empty rows with given sides */
      for( i = 0; i < nrows; ++i )
         rows.add(lhs[i], emptyVector, rhs[i]);
      spx->addRowsReal(rows);

      /* create column vectors with coefficients and bounds */
      SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                /*colnames*/,       /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);
   assert(nnonz >= 0);
   assert(ncols >= 0);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

#ifndef NDEBUG
   if ( nnonz > 0 )
   {
      /* perform check that no new rows are added - this is likely to be a mistake */
      int nrows = lpi->spx->numRowsReal();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( 0 <= ind[j] && ind[j] < nrows );
         assert( val[j] != 0.0 );
      }
   }
#endif

   SPxSCIP* spx = lpi->spx;
   try
   {
      LPColSet cols(ncols);
      DSVector colVector(ncols);
      int start;
      int last;
      int i;

      /* create column vectors with coefficients and bounds */
      for( i = 0; i < ncols; ++i )
      {
         colVector.clear();
         if( nnonz > 0 )
         {
            start = beg[i];
            last = (i == ncols-1 ? nnonz : beg[i+1]);
            colVector.add( last-start, &ind[start], &val[start] );
         }
         cols.add(obj[i], lb[i], colVector, ub[i]);
      }
      spx->addColsReal(cols);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsReal());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeColRangeReal(firstcol, lastcol) );

   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int ncols;
   int i;

   SCIPdebugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(dstat != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   ncols = lpi->spx->numColsReal();

   /* SoPlex removeCols() method deletes the columns with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < ncols; ++i )
      dstat[i] *= -1;

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeColsReal(dstat) );

   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                /*rownames*/,       /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

#ifndef NDEBUG
   if ( nnonz > 0 )
   {
      /* perform check that no new columns are added - this is likely to be a mistake */
      int ncols = lpi->spx->numColsReal();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( val[j] != 0.0 );
         assert( 0 <= ind[j] && ind[j] < ncols );
      }
   }
#endif

   try
   {
      SPxSCIP* spx = lpi->spx;
      LPRowSet rows(nrows);
      DSVector rowVector;
      int start;
      int last;
      int i;

      /* create row vectors with given sides */
      for( i = 0; i < nrows; ++i )
      {
         rowVector.clear();
         if( nnonz > 0 )
         {
            start = beg[i];
            last = (i == nrows-1 ? nnonz : beg[i+1]);
            rowVector.add( last-start, &ind[start], &val[start] );
         }
         rows.add(lhs[i], rowVector, rhs[i]);
      }
      spx->addRowsReal(rows);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsReal());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeRowRangeReal(firstrow, lastrow) );

   return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int nrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiDelRowset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   nrows = lpi->spx->numRowsReal();

   /* SoPlex removeRows() method deletes the rows with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < nrows; ++i )
      dstat[i] *= -1;

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeRowsReal(dstat) );

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiClear()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->clearLPReal() );

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices or NULL if ncols is zero */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   const SCIP_Real*      ub                  /**< values for the new upper bounds or NULL if ncols is zero */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ncols == 0 || (ind != NULL && lb != NULL && ub != NULL));
   if( ncols <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      for( i = 0; i < ncols; ++i )
      {
         assert(0 <= ind[i] && ind[i] < lpi->spx->numColsReal());

         if ( SCIPlpiIsInfinity(lpi, lb[i]) )
         {
            SCIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", ind[i]);
            return SCIP_LPERROR;
         }
         if ( SCIPlpiIsInfinity(lpi, -ub[i]) )
         {
            SCIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", ind[i]);
            return SCIP_LPERROR;
         }

         lpi->spx->changeBoundsReal(ind[i], lb[i], ub[i]);
         assert(lpi->spx->lowerReal(ind[i]) <= lpi->spx->upperReal(ind[i]) + lpi->spx->realParam(SoPlex::EPSILON_ZERO));
      }
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   if( nrows <= 0 )
      return SCIP_OKAY;

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      for( i = 0; i < nrows; ++i )
      {
         assert(0 <= ind[i] && ind[i] < lpi->spx->numRowsReal());
         lpi->spx->changeRangeReal(ind[i], lhs[i], rhs[i]);
         assert(lpi->spx->lhsReal(ind[i]) <= lpi->spx->rhsReal(ind[i]) + lpi->spx->realParam(SoPlex::EPSILON_ZERO));
      }
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgCoef()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= row && row < lpi->spx->numRowsReal());
   assert(0 <= col && col < lpi->spx->numColsReal());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->changeElementReal(row, col, newval) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, (void) lpi->spx->setIntParam(SoPlex::OBJSENSE, objsen == SCIP_OBJSEN_MINIMIZE ? SoPlex::OBJSENSE_MINIMIZE : SoPlex::OBJSENSE_MAXIMIZE ) );

   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   const int*            ind,                /**< column indices to change objective value for */
   const SCIP_Real*      obj                 /**< new objective values for columns */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      for( i = 0; i < ncols; ++i )
      {
         assert(0 <= ind[i] && ind[i] < lpi->spx->numColsReal());
         lpi->spx->changeObjReal(ind[i], obj[i]);
      }
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIPdebugMessage("calling SCIPlpiScaleRow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   try
   {
      invalidateSolution(lpi);

      assert( lpi->spx->preStrongbranchingBasisFreed() );

      /* get the row vector and the row's sides */
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
      SVector rowvec = lpi->spx->rowVectorRealInternal(row);
#else
      SVector rowvec = lpi->spx->rowVectorReal(row);
#endif
      lhs = lpi->spx->lhsReal(row);
      rhs = lpi->spx->rhsReal(row);

      /* scale the row vector */
      rowvec *= scaleval;

      /* adjust the sides */
      if( lhs > -lpi->spx->realParam(SoPlex::INFTY) )
         lhs *= scaleval;
      else if( scaleval < 0.0 )
         lhs = lpi->spx->realParam(SoPlex::INFTY);
      if( rhs < lpi->spx->realParam(SoPlex::INFTY) )
         rhs *= scaleval;
      else if( scaleval < 0.0 )
         rhs = -lpi->spx->realParam(SoPlex::INFTY);
      if( scaleval < 0.0 )
      {
         SCIP_Real oldlhs = lhs;
         lhs = rhs;
         rhs = oldlhs;
      }

      /* create the new row */
      LPRow lprow(lhs, rowvec, rhs);

      /* change the row in the LP */
      lpi->spx->changeRowReal(row, lprow);
      assert(lpi->spx->lhsReal(row) <= lpi->spx->rhsReal(row));
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   SCIPdebugMessage("calling SCIPlpiScaleCol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   try
   {
      invalidateSolution(lpi);

      assert( lpi->spx->preStrongbranchingBasisFreed() );

      /* get the col vector and the col's bounds and objective value */
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
      SVector colvec = lpi->spx->colVectorRealInternal(col);
#else
      SVector colvec = lpi->spx->colVectorReal(col);
#endif
      obj = lpi->spx->objReal(col);
      lb = lpi->spx->lowerReal(col);
      ub = lpi->spx->upperReal(col);

      /* scale the col vector */
      colvec *= scaleval;

      /* scale the objective value */
      obj *= scaleval;

      /* adjust the bounds */
      if( lb > -lpi->spx->realParam(SoPlex::INFTY) )
         lb /= scaleval;
      else if( scaleval < 0.0 )
         lb = lpi->spx->realParam(SoPlex::INFTY);
      if( ub < lpi->spx->realParam(SoPlex::INFTY) )
         ub /= scaleval;
      else if( scaleval < 0.0 )
         ub = -lpi->spx->realParam(SoPlex::INFTY);
      if( scaleval < 0.0 )
      {
         SCIP_Real oldlb = lb;
         lb = ub;
         ub = oldlb;
      }

      /* create the new col (in LPCol's constructor, the upper bound is given first!) */
      LPCol lpcol(obj, colvec, ub, lb);

      /* change the col in the LP */
      lpi->spx->changeColReal(col, lpcol);
      assert(lpi->spx->lowerReal(col) <= lpi->spx->upperReal(col));
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nrows != NULL);

   *nrows = lpi->spx->numRowsReal();

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ncols != NULL);

   *ncols = lpi->spx->numColsReal();

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetNNonz()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nnonz != NULL);

   /* SoPlex has no direct method to return the number of nonzeros, so we have to count them manually */
   *nnonz = 0;
   if( lpi->spx->numRowsReal() < lpi->spx->numColsReal() )
   {
      for( i = 0; i < lpi->spx->numRowsReal(); ++i )
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
         (*nnonz) += lpi->spx->rowVectorRealInternal(i).size();
#else
         (*nnonz) += lpi->spx->rowVectorReal(i).size();
#endif
   }
   else
   {
      for( i = 0; i < lpi->spx->numColsReal(); ++i )
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
         (*nnonz) += lpi->spx->colVectorRealInternal(i).size();
#else
         (*nnonz) += lpi->spx->colVectorReal(i).size();
#endif
   }

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiGetCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsReal());
   assert((lb != NULL && ub != NULL) || (lb == NULL && ub == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   if( lb != NULL )
   {
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
      if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) )
      {
         DVector lbvec(lpi->spx->numColsReal());
         DVector ubvec(lpi->spx->numColsReal());
         lpi->spx->getLowerReal(lbvec);
         lpi->spx->getUpperReal(ubvec);
         for( i = firstcol; i <= lastcol; ++i )
         {
            lb[i-firstcol] = lbvec[i];
            ub[i-firstcol] = ubvec[i];
         }
      }
      else
      {
         const Vector& lbvec = lpi->spx->lowerRealInternal();
         const Vector& ubvec = lpi->spx->upperRealInternal();
         for( i = firstcol; i <= lastcol; ++i )
         {
            lb[i-firstcol] = lbvec[i];
            ub[i-firstcol] = ubvec[i];
         }
      }
#else
      const Vector& lbvec = lpi->spx->lowerReal();
      const Vector& ubvec = lpi->spx->upperReal();

      for( i = firstcol; i <= lastcol; ++i )
      {
         lb[i-firstcol] = lbvec[i];
         ub[i-firstcol] = ubvec[i];
      }
#endif
   }

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstcol; i <= lastcol; ++i )
      {
         beg[i-firstcol] = *nnonz;

#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
         if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) )
         {
            DSVector cvec;
            lpi->spx->getColVectorReal(i, cvec);
            for( j = 0; j < cvec.size(); ++j )
            {
               ind[*nnonz] = cvec.index(j);
               val[*nnonz] = cvec.value(j);
               (*nnonz)++;
            }
         }
         else
         {
            const SVector& cvec = lpi->spx->colVectorRealInternal(i);
            for( j = 0; j < cvec.size(); ++j )
            {
               ind[*nnonz] = cvec.index(j);
               val[*nnonz] = cvec.value(j);
               (*nnonz)++;
            }
         }
#else
         const SVector& cvec = lpi->spx->colVectorReal(i);
         for( j = 0; j < cvec.size(); ++j )
         {
            ind[*nnonz] = cvec.index(j);
            val[*nnonz] = cvec.value(j);
            (*nnonz)++;
         }
#endif
      }
   }

   return SCIP_OKAY;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiGetRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsReal());
   assert((lhs != NULL && rhs != NULL) || (lhs == NULL && rhs == NULL));
   assert((nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL));

   if( lhs != NULL )
   {
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
      if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) )
      {
         DVector lhsvec(lpi->spx->numRowsReal());
         DVector rhsvec(lpi->spx->numRowsReal());
         lpi->spx->getLhsReal(lhsvec);
         lpi->spx->getRhsReal(rhsvec);
         for( i = firstrow; i <= lastrow; ++i )
         {
            lhs[i-firstrow] = lhsvec[i];
            rhs[i-firstrow] = rhsvec[i];
         }
      }
      else
      {
         const Vector& lhsvec = lpi->spx->lhsRealInternal();
         const Vector& rhsvec = lpi->spx->rhsRealInternal();
         for( i = firstrow; i <= lastrow; ++i )
         {
            lhs[i-firstrow] = lhsvec[i];
            rhs[i-firstrow] = rhsvec[i];
         }
      }
#else
      const Vector& lhsvec = lpi->spx->lhsReal();
      const Vector& rhsvec = lpi->spx->rhsReal();
      for( i = firstrow; i <= lastrow; ++i )
      {
         lhs[i-firstrow] = lhsvec[i];
         rhs[i-firstrow] = rhsvec[i];
      }
#endif
   }

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstrow; i <= lastrow; ++i )
      {
         beg[i-firstrow] = *nnonz;

#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
         if( lpi->spx->boolParam(SoPlex::PERSISTENTSCALING) )
         {
            DSVector rvec;
            lpi->spx->getRowVectorReal(i, rvec);
            for( j = 0; j < rvec.size(); ++j )
            {
               ind[*nnonz] = rvec.index(j);
               val[*nnonz] = rvec.value(j);
               (*nnonz)++;
            }
         }
         else
         {
            const SVector& rvec = lpi->spx->rowVectorRealInternal(i);
            for( j = 0; j < rvec.size(); ++j )
            {
               ind[*nnonz] = rvec.index(j);
               val[*nnonz] = rvec.value(j);
               (*nnonz)++;
            }
         }
#else
         const SVector& rvec = lpi->spx->rowVectorReal(i);
         for( j = 0; j < rvec.size(); ++j )
         {
            ind[*nnonz] = rvec.index(j);
            val[*nnonz] = rvec.value(j);
            (*nnonz)++;
         }
#endif
      }
   }

   return SCIP_OKAY;
}

/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( colnames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsReal() );

   SCIPdebugMessage("getting column names %d to %d\n", firstcol, lastcol);

//    lpi->spx->getColNames(firstcol, lastcol, colnames, namestorage, namestoragesize, storageleft);

   return SCIP_OKAY;
}

/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( rownames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsReal() );

   SCIPdebugMessage("getting row names %d to %d\n", firstrow, lastrow);

//    lpi->spx->getRowNames(firstrow, lastrow, rownames, namestorage, namestoragesize, storageleft);

   return SCIP_OKAY;
}

/** gets objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(objsen != NULL);

   *objsen = (lpi->spx->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE) ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE;

   return SCIP_OKAY;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsReal());
   assert(vals != NULL);

   for( i = firstcol; i <= lastcol; ++i )
      vals[i-firstcol] = lpi->spx->objReal(i);

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->numColsReal());

   for( i = firstcol; i <= lastcol; ++i )
   {
      if( lbs != NULL )
         lbs[i-firstcol] = lpi->spx->lowerReal(i);
      if( ubs != NULL )
         ubs[i-firstcol] = lpi->spx->upperReal(i);
   }

   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->numRowsReal());

   for( i = firstrow; i <= lastrow; ++i )
   {
      if( lhss != NULL )
         lhss[i-firstrow] = lpi->spx->lhsReal(i);
      if( rhss != NULL )
         rhss[i-firstrow] = lpi->spx->rhsReal(i);
   }

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetCoef()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= col && col < lpi->spx->numColsReal());
   assert(0 <= row && row < lpi->spx->numRowsReal());
   assert(val != NULL);

#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
   *val = lpi->spx->coefReal(row, col);
#else
   *val = lpi->spx->colVectorReal(col)[row];
#endif

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves LP -- used for both, primal and dual simplex, because SoPlex doesn't distinct the two cases */
static
SCIP_RETCODE spxSolve(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );

   SPxOut::Verbosity verbosity;
   /* store and set verbosity */
   verbosity = lpi->spx->spxout.getVerbosity();
   lpi->spx->spxout.setVerbosity((SPxOut::Verbosity)(lpi->spx->getLpInfo() ? SOPLEX_VERBLEVEL : 0));

   SCIPdebugMessage("calling SoPlex solve(): %d cols, %d rows\n", lpi->spx->numColsReal(), lpi->spx->numRowsReal());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

#ifdef WITH_LPSCHECK
   lpi->spx->setDoubleCheck(CHECK_SPXSOLVE);
#endif

   /* delete starting basis if solving from scratch */
   if( lpi->spx->getFromScratch() )
   {
      try
      {
         lpi->spx->clearBasis();
      }
#ifndef NDEBUG
      catch(const SPxException& x)
      {
         std::string s = x.what();
         SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
      catch(const SPxException&)
      {
#endif
         assert( lpi->spx->status() != SPxSolver::OPTIMAL );
         return SCIP_LPERROR;
      }
   }
   assert(!lpi->spx->getFromScratch() || lpi->spx->status() == SPxSolver::NO_PROBLEM);

   SPxSolver::Status status = lpi->spx->doSolve();
   SCIPdebugMessage(" -> SoPlex status: %d, basis status: %d\n", lpi->spx->status(), lpi->spx->basisStatus());
   lpi->solved = TRUE;

   /* restore verbosity */
   lpi->spx->spxout.setVerbosity(verbosity);

   switch( status )
   {
   case SPxSolver::ABORT_TIME:
   case SPxSolver::ABORT_ITER:
   case SPxSolver::ABORT_VALUE:
   case SPxSolver::SINGULAR:
   case SPxSolver::REGULAR:
   case SPxSolver::UNKNOWN:
   case SPxSolver::OPTIMAL:
   case SPxSolver::UNBOUNDED:
   case SPxSolver::INFEASIBLE:
      return SCIP_OKAY;
   default:
      return SCIP_LPERROR;
   }  /*lint !e788*/
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolvePrimal()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   (void) lpi->spx->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_PRIMAL);
   return spxSolve(lpi);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolveDual()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   (void) lpi->spx->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_DUAL);
   return spxSolve(lpi);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   SCIPdebugMessage("calling SCIPlpiSolveBarrier()\n");

   /* Since SoPlex does not support barrier we switch to DUAL */
   return SCIPlpiSolveDual(lpi);
}

/** start strong branching - call before any strongbranching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   lpi->spx->savePreStrongbranchingBasis();

   return SCIP_OKAY;
}

/** end strong branching - call after any strongbranching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( ! lpi->spx->preStrongbranchingBasisFreed() );
   lpi->spx->restorePreStrongbranchingBasis();
   lpi->spx->freePreStrongbranchingBasis();

   return SCIP_OKAY;
}

/** performs strong branching iterations on one arbitrary candidate */
static
SCIP_RETCODE lpiStrongbranch(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   SPxSCIP* spx;
   SPxSolver::Status status;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   bool fromparentbasis;
   bool error;
   int oldItlim;
   SPxOut::Verbosity verbosity;

   /* store and set verbosity */
   verbosity = lpi->spx->spxout.getVerbosity();
   lpi->spx->spxout.setVerbosity((SPxOut::Verbosity)(lpi->spx->getLpInfo() ? SOPLEX_VERBLEVEL : 0));

   SCIPdebugMessage("calling SCIPlpiStrongbranch() on variable %d (%d iterations)\n", col, itlim);

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   /*  assert(down != NULL);
    * assert(up != NULL); temporary hack for cloud branching */
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   spx = lpi->spx;
#ifndef STRONGBRANCH_RESTOREBASIS
   fromparentbasis = false;
#endif
   error = false;
   oldItlim = spx->intParam(SoPlex::ITERLIMIT);

   /* get current bounds of column */
   oldlb = spx->lowerReal(col);
   oldub = spx->upperReal(col);

   *downvalid = FALSE;
   *upvalid = FALSE;

   if( iter != NULL )
      *iter = 0;

   /* set the algorithm type to use dual simplex */
   (void) spx->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_DUAL);

   /* down branch */
   newub = EPSCEIL(psol-1.0, lpi->spx->feastol());
   if( newub >= oldlb - 0.5 && down != NULL )
   {
      SCIPdebugMessage("strong branching down on x%d (%g) with %d iterations\n", col, psol, itlim);

      spx->changeUpperReal(col, newub);
      assert(spx->lowerReal(col) <= spx->upperReal(col));

      (void) spx->setIntParam(SoPlex::ITERLIMIT, itlim);
      do
      {
#ifdef WITH_LPSCHECK
         spx->setDoubleCheck(CHECK_SPXSTRONGBRANCH);
#endif
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
         status =  spx->optimize();
#else
         status = spx->solve();
#endif
         SCIPdebugMessage(" --> Terminate with status %d\n", status);
         switch( status )
         {
         case SPxSolver::OPTIMAL:
            *down = spx->objValueReal();
            *downvalid = TRUE;
            SCIPdebugMessage(" --> Terminate with value %f\n", *down);
            break;
         case SPxSolver::ABORT_TIME: /* SoPlex does not return a proven dual bound, if it is aborted */
         case SPxSolver::ABORT_ITER:
         case SPxSolver::ABORT_CYCLING:
            *down = spx->objValueReal();
            break;
         case SPxSolver::ABORT_VALUE:
         case SPxSolver::INFEASIBLE:
            *down = spx->getObjLimit();
            *downvalid = TRUE;
            break;
         default:
            error = true;
            break;
         }  /*lint !e788*/
         if( iter != NULL )
            (*iter) += spx->numIterations();

#ifdef STRONGBRANCH_RESTOREBASIS
         /* we restore the pre-strong-branching basis by default (and don't solve again) */
         assert( ! spx->preStrongbranchingBasisFreed() );
         spx->restorePreStrongbranchingBasis();
         fromparentbasis = false;
#else
         /* if cycling or singular basis occured and we started not from the pre-strong-branching basis, then we restore the
          * pre-strong-branching basis and try again with reduced iteration limit */
         if( (status == SPxSolver::ABORT_CYCLING || status == SPxSolver::SINGULAR) && !fromparentbasis && spx->numIterations() < itlim )
         {
            SCIPdebugMessage(" --> Repeat strong branching down with %d iterations after restoring basis\n",
                             itlim - spx->numIterations());
            spx->setIntParam(SoPlex::ITERLIMIT, itlim - spx->numIterations());
            assert( ! spx->hasPreStrongbranchingBasis() );
            spx->restorePreStrongbranchingBasis();
            fromparentbasis = true;
            error = false;
         }
         /* otherwise don't solve again */
         else
            fromparentbasis = false;
#endif
      }
      while( fromparentbasis );

      spx->changeUpperReal(col, oldub);
      assert(spx->lowerReal(col) <= spx->upperReal(col));
   }
   else if( down != NULL )
   {
      *down = spx->getObjLimit();
      *downvalid = TRUE;
   }
   else
      *downvalid = TRUE;

   /* up branch */
   if( !error )
   {
      newlb = EPSFLOOR(psol+1.0, lpi->spx->feastol());
      if( newlb <= oldub + 0.5 && up != NULL )
      {
         SCIPdebugMessage("strong branching  up  on x%d (%g) with %d iterations\n", col, psol, itlim);

         spx->changeLowerReal(col, newlb);
         assert(spx->lowerReal(col) <= spx->upperReal(col));

         (void) spx->setIntParam(SoPlex::ITERLIMIT, itlim);
         do
         {
#ifdef WITH_LPSCHECK
            spx->setDoubleCheck(CHECK_SPXSTRONGBRANCH);
#endif
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
            status = spx->optimize();
#else
            status = spx->solve();
#endif
            SCIPdebugMessage(" --> Terminate with status %d\n", status);
            switch( status )
            {
            case SPxSolver::OPTIMAL:
               *up = spx->objValueReal();
               *upvalid = TRUE;
               SCIPdebugMessage(" --> Terminate with value %f\n", spx->objValueReal());
               break;
            case SPxSolver::ABORT_TIME: /* SoPlex does not return a proven dual bound, if it is aborted */
            case SPxSolver::ABORT_ITER:
            case SPxSolver::ABORT_CYCLING:
               *up = spx->objValueReal();
               break;
            case SPxSolver::ABORT_VALUE:
            case SPxSolver::INFEASIBLE:
               *up = spx->getObjLimit();
               *upvalid = TRUE;
               break;
            default:
               error = true;
               break;
            }  /*lint !e788*/
            if( iter != NULL )
               (*iter) += spx->numIterations();

#ifdef STRONGBRANCH_RESTOREBASIS
            /* we restore the pre-strong-branching basis by default (and don't solve again) */
            assert( ! spx->preStrongbranchingBasisFreed() );
            spx->restorePreStrongbranchingBasis();
            fromparentbasis = false;
#else
            /* if cycling or singular basis occured and we started not from the pre-strong-branching basis, then we restore the
             * pre-strong-branching basis and try again with reduced iteration limit */
            else if( (status == SPxSolver::ABORT_CYCLING || status == SPxSolver::SINGULAR) && !fromparentbasis && spx->numIterations() < itlim )
            {
               SCIPdebugMessage(" --> Repeat strong branching  up  with %d iterations after restoring basis\n", itlim - spx->numIterations());
               assert( ! spx->hasPreStrongbranchingBasis() );
               spx->restorePreStrongbranchingBasis();
               spx->setIntParam(SoPlex::ITERLIMIT, itlim - spx->numIterations());
               error = false;
               fromparentbasis = true;
            }
            /* otherwise don't solve again */
            else
               fromparentbasis = false;
#endif
         }
         while( fromparentbasis );

         spx->changeLowerReal(col, oldlb);
         assert(spx->lowerReal(col) <= spx->upperReal(col));
      }
      else if( up != NULL )
      {
         *up = spx->getObjLimit();
         *upvalid = TRUE;
      }
      else
         *upvalid = TRUE;
   }

   /* reset old iteration limit */
   (void) spx->setIntParam(SoPlex::ITERLIMIT, oldItlim);

   /* restore verbosity */
   lpi->spx->spxout.setVerbosity(verbosity);

   if( error )
   {
      SCIPdebugMessage("SCIPlpiStrongbranch() returned SoPlex status %d\n", int(status));  /*lint !e644*/
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   /* pass call on to lpiStrongbranch() */
   retcode = lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter);

   /* pass SCIP_LPERROR to SCIP without a back trace */
   if( retcode == SCIP_LPERROR )
      return SCIP_LPERROR;

   /* evaluate retcode */
   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (int j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      retcode = lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter);

      /* pass SCIP_LPERROR to SCIP without a back trace */
      if( retcode == SCIP_LPERROR )
         return SCIP_LPERROR;

      /* evaluate retcode */
      SCIP_CALL( retcode );
   }
   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   /* pass call on to lpiStrongbranch() */
   retcode = lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter);

   /* pass SCIP_LPERROR to SCIP without a back trace */
   if( retcode == SCIP_LPERROR )
      return SCIP_LPERROR;

   /* evaluate retcode */
   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (int j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      retcode = lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter);

      /* pass SCIP_LPERROR to SCIP without a back trace */
      if( retcode == SCIP_LPERROR )
         return SCIP_LPERROR;

      /* evaluate retcode */
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   return lpi->solved;
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolFeasibility()\n");

   assert(lpi != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   *primalfeasible = SCIPlpiIsPrimalFeasible(lpi);
   *dualfeasible = SCIPlpiIsDualFeasible(lpi);

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return lpi->spx->hasPrimalRay();
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert(lpi->spx->status() != SPxSolver::UNBOUNDED || lpi->spx->basisStatus() == SPxBasis::UNBOUNDED);

   /* if SoPlex returns unbounded, this may only mean that an unbounded ray is available, not necessarily a primal
    * feasible point; hence we have to check the perturbation
    */
   return lpi->spx->status() == SPxSolver::UNBOUNDED;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SPxBasis::SPxStatus basestatus;

   SCIPdebugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   basestatus = lpi->spx->basisStatus();

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   assert(basestatus == SPxBasis::OPTIMAL || lpi->spx->status() != SPxSolver::OPTIMAL);

   return basestatus == SPxBasis::OPTIMAL || basestatus == SPxBasis::PRIMAL;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return lpi->spx->hasDualFarkas();
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return lpi->spx->status() == SPxSolver::INFEASIBLE && lpi->spx->basisStatus() == SPxBasis::DUAL;
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   assert(lpi->spx->basisStatus() == SPxBasis::OPTIMAL || lpi->spx->status() != SPxSolver::OPTIMAL);

   return (lpi->spx->basisStatus() == SPxBasis::OPTIMAL) || lpi->spx->basisStatus() == SPxBasis::DUAL;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert((lpi->spx->basisStatus() == SPxBasis::OPTIMAL)
      == (SCIPlpiIsPrimalFeasible(lpi) && SCIPlpiIsDualFeasible(lpi)));

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   return (lpi->spx->basisStatus() == SPxBasis::OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( lpi->spx->status() == SPxSolver::ERROR || lpi->spx->status() == SPxSolver::SINGULAR )
      return FALSE;

   /* only if we have a regular basis and the condition limit is set, we compute the condition number of the basis;
    * everything above the specified threshold is then counted as instable
    */
   if( lpi->checkcondition && (SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi)) )
   {
      SCIP_RETCODE retcode;
      SCIP_Real kappa;

      retcode = SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_ESTIMCONDITION, &kappa);
      if( retcode != SCIP_OKAY )
      {
         SCIPABORT();
      }
      assert(kappa != SCIP_INVALID); /*lint !e777*/

      if( kappa > lpi->conditionlimit )
         return FALSE;
   }

   return TRUE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::ABORT_VALUE);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::ABORT_ITER);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->status() == SPxSolver::ABORT_TIME);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return static_cast<int>(lpi->spx->status());
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiIgnoreInstability()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(success != NULL);

   /* instable situations cannot be ignored */
   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjval()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(objval != NULL);

   *objval = lpi->spx->objValueReal();

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( objval != NULL )
      *objval = lpi->spx->objValueReal();

   try
   {
      if( primsol != NULL )
      {
         Vector tmp(lpi->spx->numColsReal(), primsol);
         (void)lpi->spx->getPrimalReal(tmp);
      }
      if( dualsol != NULL )
      {
         Vector tmp(lpi->spx->numRowsReal(), dualsol);
         (void)lpi->spx->getDualReal(tmp);
      }
      if( activity != NULL )
      {
         Vector tmp(lpi->spx->numRowsReal(), activity);
         (void)lpi->spx->getSlacksReal(tmp);  /* in SoPlex, the activities are called "slacks" */
      }
      if( redcost != NULL )
      {
         Vector tmp(lpi->spx->numColsReal(), redcost);
         (void)lpi->spx->getRedCostReal(tmp);
      }
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiGetPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpi->spx->hasPrimalRay());
   assert(ray != NULL);

   try
   {
      Vector tmp(lpi->spx->numColsReal(), ray);
      (void)lpi->spx->getPrimalRayReal(tmp);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpi->spx->hasDualFarkas());
   assert(dualfarkas != NULL);

   try
   {
      Vector tmp(lpi->spx->numRowsReal(), dualfarkas);
      (void)lpi->spx->getDualFarkasReal(tmp);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetIterations()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(iterations != NULL);

   *iterations = lpi->spx->numIterations();

   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for @p quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRealSolQuality()\n");

   assert(lpi != NULL);
   assert(quality != NULL);

   bool success;

   SCIPdebugMessage("requesting solution quality from SoPlex: quality %d\n", qualityindicator);

   switch( qualityindicator )
   {
      case SCIP_LPSOLQUALITY_ESTIMCONDITION:
         success = lpi->spx->getEstimatedCondition(*quality);
         break;

      case SCIP_LPSOLQUALITY_EXACTCONDITION:
         success = lpi->spx->getExactCondition(*quality);
         break;

      default:
         SCIPerrorMessage("Solution quality %d unknown.\n", qualityindicator);
         return SCIP_INVALIDDATA;
   }

   if( !success )
   {
      SCIPdebugMessage("problem computing condition number\n");
      *quality = SCIP_INVALID;
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */


/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   if( rstat != NULL )
   {
      for( i = 0; i < lpi->spx->numRowsReal(); ++i )
      {
         switch( lpi->spx->basisRowStatus(i) )
         {
         case SPxSolver::BASIC:
            rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
            break;
         case SPxSolver::FIXED:
         case SPxSolver::ON_LOWER:
            rstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_UPPER:
            rstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
            break;
         case SPxSolver::ZERO:
            SCIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
            return SCIP_LPERROR;
         case SPxSolver::UNDEFINED:
         default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   if( cstat != NULL )
   {
      for( i = 0; i < lpi->spx->numColsReal(); ++i )
      {
//         SCIP_Real val = 0.0;
         switch( lpi->spx->basisColStatus(i) )
         {
         case SPxSolver::BASIC:
            cstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
            break;
         case SPxSolver::FIXED:
         /* Get reduced cost estimation. If the estimation is not correct this should not hurt:
         * If the basis is loaded into SoPlex again, the status is converted to FIXED again; in
         * this case there is no problem at all. If the basis is saved and/or used in some other
         * solver, it usually is very cheap to perform the pivots necessary to get an optimal
         * basis.
         * @todo implement getRedCostEst()
         * */
//          SCIP_CALL( getRedCostEst(lpi->spx, i, &val) );
//            if( val < 0.0 )  /* reduced costs < 0 => UPPER  else => LOWER */
//               cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
//            else
            cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_LOWER:
            cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_UPPER:
            cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
            break;
         case SPxSolver::ZERO:
            cstat[i] = SCIP_BASESTAT_ZERO; /*lint !e641*/
            break;
         case SPxSolver::UNDEFINED:
         default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const int*            cstat,              /**< array with column basis status */
   const int*            rstat               /**< array with row basis status */
   )
{
   int i;
   int ncols;
   int nrows;

   SCIPdebugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );

   assert(cstat != NULL || ncols == 0);
   assert(rstat != NULL || nrows == 0);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   invalidateSolution(lpi);

   DataArray<SPxSolver::VarStatus>& _colstat = lpi->spx->colStat();
   DataArray<SPxSolver::VarStatus>& _rowstat = lpi->spx->rowStat();

   _colstat.reSize(ncols);
   _rowstat.reSize(nrows);

   for( i = 0; i < nrows; ++i )
   {
      switch( rstat[i] ) /*lint !e613*/
      {
      case SCIP_BASESTAT_LOWER:
         _rowstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         _rowstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         _rowstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         SCIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
         return SCIP_LPERROR; /*lint !e429*/
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   for( i = 0; i < ncols; ++i )
   {
      switch( cstat[i] ) /*lint !e613*/
      {
      case SCIP_BASESTAT_LOWER:
         _colstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         _colstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         _colstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         _colstat[i] = SPxSolver::ZERO;
         break;
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->setBasis(_rowstat.get_ptr(), _colstat.get_ptr()) );
   lpi->spx->freePreStrongbranchingBasis();

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBasisInd()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(bind != NULL);

   assert(lpi->spx->preStrongbranchingBasisFreed());

   lpi->spx->getBasisInd(bind);

   return SCIP_OKAY;
}


/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvRow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpi->spx->preStrongbranchingBasisFreed());
   assert(coef != NULL);

   assert(r >= 0);
   assert(r < lpi->spx->numRowsReal());

   if( ! lpi->spx->getBasisInverseRowReal(r, coef, inds, ninds) )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvCol()\n");

   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( lpi->spx->preStrongbranchingBasisFreed() );
   assert(coef != NULL);

   if( ! lpi->spx->getBasisInverseColReal(c, coef, inds, ninds) )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< vector to return coefficients */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   SCIP_Real* buf;
   SCIP_Real* binv;
   int nrows;
   int ncols;
   int c;

   SCIPdebugMessage("calling SCIPlpiGetBInvARow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert( lpi->spx->preStrongbranchingBasisFreed() );
   assert(coef != NULL);

   nrows = lpi->spx->numRowsReal();
   ncols = lpi->spx->numColsReal();
   buf = NULL;

   /* get (or calculate) the row in B^-1 */
   if( binvrow == NULL )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&buf, nrows) );
      SCIP_CALL( SCIPlpiGetBInvRow(lpi, r, buf, inds, ninds) );
      binv = buf;
   }
   else
      binv = const_cast<SCIP_Real*>(binvrow);

   assert(binv != NULL);

   /* mark sparsity pattern as invalid */
   if( ninds != NULL )
      *ninds = -1;

   // @todo exploit sparsity in binv by looping over nrows
   /* calculate the scalar product of the row in B^-1 and A */
   Vector binvvec(nrows, binv);

#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
   /* temporary unscaled column of A */
   DSVector acol;
#endif

   for( c = 0; c < ncols; ++c )
   {
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
      lpi->spx->getColVectorReal(c, acol);
      coef[c] = binvvec * acol;  /* scalar product */ /*lint !e1702*/
#else
      coef[c] = binvvec * lpi->spx->colVectorReal(c);  /* scalar product */ /*lint !e1702*/
#endif
   }

   /* free memory if it was temporarily allocated */
   BMSfreeMemoryArrayNull(&buf);

   return SCIP_OKAY;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef,               /**< vector to return coefficients */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{  /*lint --e{715}*/
   /* create a new uninitialized full vector */
   DVector col(lpi->spx->numRowsReal());

#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
   /* temporary sparse vector used for unscaling (memory is automatically enlarged) */
   DSVector colsparse;
#endif

   SCIPdebugMessage("calling SCIPlpiGetBInvACol()\n");

   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( lpi->spx->preStrongbranchingBasisFreed() );
   assert(coef != NULL);

   /* extract column c of A */
   assert(c >= 0);
   assert(c < lpi->spx->numColsReal());

   /* @todo implement this with sparse vectors */
   /* mark sparsity pattern as invalid */
   if( ninds != NULL )
      *ninds = -1;

   /* col needs to be cleared because copying colVectorReal only regards nonzeros */
   col.clear();

#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 4)
   lpi->spx->getColVectorReal(c, colsparse);
   /* the copy is necessary to transform the sparse column into a dense vector */
   col = colsparse;
#else
   col = lpi->spx->colVectorReal(c);
#endif

   /* solve */
   if( ! lpi->spx->getBasisInverseTimesVecReal(col.get_ptr(), coef) )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/**@} */




/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   int ncols;
   int nrows;

   SCIPdebugMessage("calling SCIPlpiGetState()\n");

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   ncols = lpi->spx->numColsReal();
   nrows = lpi->spx->numRowsReal();
   assert(ncols >= 0);
   assert(nrows >= 0);

   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information */
   SCIP_CALL( SCIPlpiGetBase(lpi, lpi->cstat, lpi->rstat) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->cstat, lpi->rstat);

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           /*blkmem*/,         /**< block memory */
   const SCIP_LPISTATE*  lpistate            /**< LPi state information (like basis information), or NULL */
   )
{
   int lpncols;
   int lpnrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiSetState()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);
   /* assert(blkmem != NULL); */

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   lpncols = lpi->spx->numColsReal();
   lpnrows = lpi->spx->numRowsReal();
   assert(lpistate->ncols <= lpncols);
   assert(lpistate->nrows <= lpnrows);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, lpncols) );
   SCIP_CALL( ensureRstatMem(lpi, lpnrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < lpncols; ++i )
   {
      SCIP_Real bnd = lpi->spx->lowerReal(i);
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         bnd = lpi->spx->lowerReal(i);
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            /* variable is free */
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /*lint !e641*/
         else
            /* use finite upper bound */
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
      }
      else
         /* use finite lower bound */
         lpi->cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
   }
   for( i = lpistate->nrows; i < lpnrows; ++i )
      lpi->rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/

   /* load basis information */
   SCIP_CALL( SCIPlpiSetBase(lpi, lpi->cstat, lpi->rstat) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiClearState()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   try
   {
      lpi->spx->clearBasis();
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      assert( lpi->spx->status() != SPxSolver::OPTIMAL );
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiFreeState()\n");

   assert(lpi != NULL);
   assert(lpistate != NULL);
   assert(blkmem != NULL);

   if ( *lpistate != NULL )
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL */
   )
{  /*lint --e{715}*/
   assert(lpi != NULL);
   return TRUE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,               /**< LP interface structure */
   const char*           fname              /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadState()\n");
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   bool success;
   SOPLEX_TRY( lpi->messagehdlr, success = lpi->spx->readBasisFile(fname, 0, 0) );

   return success ? SCIP_OKAY : SCIP_LPERROR;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,            /**< LP interface structure */
   const char*           fname           /**< file name */
   )
{
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);
   SCIPdebugMessage("calling SCIPlpiWriteState()\n");

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   bool res;
   SOPLEX_TRY( lpi->messagehdlr, res = lpi->spx->writeBasisFile(fname, 0, 0) );

   if ( ! res )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Pricing Norms Methods
 */

/**@name LP Pricing Norms Methods */
/**@{ */

/** stores LPi pricing norms information
 *  @todo should we store norm information?
 */
SCIP_RETCODE SCIPlpiGetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{  /*lint --e{715}*/
#if ((SOPLEX_VERSION == 201 && SOPLEX_SUBVERSION >= 3) || SOPLEX_VERSION > 201)
   int nrows;
   int ncols;

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpinorms != NULL);

   lpi->spx->getNdualNorms(nrows, ncols);

   if( nrows == 0 && ncols == 0)
   {
      (*lpinorms = NULL);
      return SCIP_OKAY;
   }

   /* allocate lpinorms data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpinorms) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpinorms)->norms, nrows + ncols) );
   (*lpinorms)->nrows = 0;
   (*lpinorms)->ncols = 0;

   SCIPdebugMessage("storing SoPlex LPi pricing norms in %p (%d rows, %d cols)\n", (void *) *lpinorms, nrows, ncols);

   if( !lpi->spx->getDualNorms((*lpinorms)->nrows, (*lpinorms)->ncols, (*lpinorms)->norms) )
   {
      SCIPdebugMessage("freeing norms at %p\n", (void *) *lpinorms);
      BMSfreeBlockMemoryArray(blkmem, &(*lpinorms)->norms, nrows + ncols);
      BMSfreeBlockMemory(blkmem, lpinorms);
      assert(*lpinorms == NULL);
   }
#ifndef NDEBUG
   else
   {
      assert(nrows == (*lpinorms)->nrows);
      assert(ncols == (*lpinorms)->ncols);
   }
#endif
#else
   (*lpinorms) = NULL;
#endif

   return SCIP_OKAY;
}

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetNorms()
 */
SCIP_RETCODE SCIPlpiSetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPINORMS*  lpinorms            /**< LPi pricing norms information, or NULL */
   )
{  /*lint --e{715}*/
#if ((SOPLEX_VERSION == 201 && SOPLEX_SUBVERSION >= 3) || SOPLEX_VERSION > 201)
   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* if there was no pricing norms information available, the LPi norms were not stored */
   if( lpinorms == NULL )
      return SCIP_OKAY;

   assert(lpinorms->nrows <= lpi->spx->numRowsReal());
   assert(lpinorms->ncols <= lpi->spx->numColsReal());

   if( lpinorms->nrows == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("loading LPi simplex norms %p (%d rows, %d cols) into SoPlex LP with %d rows and %d cols\n",
      (const void *) lpinorms, lpinorms->nrows, lpinorms->ncols, lpi->spx->numRowsReal(), lpi->spx->numColsReal());

   (void) lpi->spx->setDualNorms(lpinorms->nrows, lpinorms->ncols, lpinorms->norms);
#endif

   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   )
{  /*lint --e{715}*/
#if ((SOPLEX_VERSION == 201 && SOPLEX_SUBVERSION >= 3) || SOPLEX_VERSION > 201)
   assert(lpi != NULL);
   assert(lpinorms != NULL);

   SCIPdebugMessage("freeing norms at %p\n", (void *) *lpinorms);

   BMSfreeBlockMemoryArray(blkmem, &(*lpinorms)->norms, (*lpinorms)->nrows + (*lpinorms)->ncols);
   BMSfreeBlockMemory(blkmem, lpinorms);
   assert(*lpinorms == NULL);
#endif

   return SCIP_OKAY;
}

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   int scaleparam;

   SCIPdebugMessage("calling SCIPlpiGetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ival != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = lpi->spx->getFromScratch();
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = lpi->spx->getLpInfo();
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = lpi->spx->intParam(SoPlex::ITERLIMIT);
      break;
   case SCIP_LPPAR_PRESOLVING:
      *ival = lpi->spx->intParam(SoPlex::SIMPLIFIER) == SoPlex::SIMPLIFIER_AUTO;
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int) lpi->pricing;
      break;
   case SCIP_LPPAR_SCALING:
      scaleparam = lpi->spx->intParam(SoPlex::SCALER);

      if( scaleparam == SoPlex::SCALER_OFF )
         *ival = 0;
      else if( scaleparam == SoPlex::SCALER_BIEQUI )
         *ival = 1;
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 2)
      else
      {
         assert(scaleparam == SoPlex::SCALER_LEASTSQ);
         *ival = 2;
      }
#else
      else
      {
         assert(scaleparam == SoPlex::SCALER_GEO8);
         *ival = 2;
      }
#endif
      break;
#if SOPLEX_VERSION >= 201
   case SCIP_LPPAR_TIMING:
      *ival = (int) (lpi->spx->intParam(SoPlex::TIMER));
      break;
#endif
#if SOPLEX_VERSION >= 230 || (SOPLEX_VERSION == 220 && SOPLEX_SUBVERSION >= 3)
   case SCIP_LPPAR_RANDOMSEED:
      *ival = (int) lpi->spx->randomSeed();
      break;
#endif
#if SOPLEX_APIVERSION >= 1
   case SCIP_LPPAR_REFACTOR:
      *ival = (int) lpi->spx->intParam(SoPlex::FACTOR_UPDATE_MAX);
      break;
#endif
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setFromScratch(bool(ival));
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setLpInfo(bool(ival));
      break;
   case SCIP_LPPAR_LPITLIM:
      assert(ival >= -1);
      (void) lpi->spx->setIntParam(SoPlex::ITERLIMIT, ival);
      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      (void) lpi->spx->setIntParam(SoPlex::SIMPLIFIER, (ival ? SoPlex::SIMPLIFIER_AUTO : SoPlex::SIMPLIFIER_OFF));
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = (SCIP_PRICING)ival;
      switch( lpi->pricing )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_AUTO:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_AUTO);
         break;
      case SCIP_PRICING_FULL:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP);
         break;
      case SCIP_PRICING_PARTIAL:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_PARMULT);
         break;
      case SCIP_PRICING_STEEP:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP);
         break;
      case SCIP_PRICING_STEEPQSTART:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_QUICKSTEEP);
         break;
      case SCIP_PRICING_DEVEX:
         (void) lpi->spx->setIntParam(SoPlex::PRICER, SoPlex::PRICER_DEVEX);
         break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival >= 0 && ival <= 2);
      if( ival == 0 )
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_OFF);
      else if( ival == 1 )
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_BIEQUI);
      else
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 2)
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_LEASTSQ);
#else
         (void) lpi->spx->setIntParam(SoPlex::SCALER, SoPlex::SCALER_GEO8);
#endif

      break;
#if SOPLEX_VERSION >= 201
   case SCIP_LPPAR_TIMING:
      assert(ival >= 0 && ival < 3);
      (void) lpi->spx->setIntParam(SoPlex::TIMER, ival);
      break;
#endif
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION == 221 && SOPLEX_SUBVERSION >= 3)
   case SCIP_LPPAR_RANDOMSEED:
      lpi->spx->setRandomSeed((unsigned long)(long)ival);
      break;
#endif
#if SOPLEX_VERSION > 221 || (SOPLEX_VERSION >= 221 && SOPLEX_SUBVERSION >= 3)
   case SCIP_LPPAR_POLISHING:
      assert(ival >= 0 && ival < 3);
      (void) lpi->spx->setIntParam(SoPlex::SOLUTION_POLISHING, ival);
      break;
#endif
#if SOPLEX_APIVERSION >= 1
   case SCIP_LPPAR_REFACTOR:
      assert(ival >= 0);
      (void) lpi->spx->setIntParam(SoPlex::FACTOR_UPDATE_MAX, ival);
      break;
#endif

   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(dval != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = lpi->spx->feastol();
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      *dval = lpi->spx->opttol();
      break;
   case SCIP_LPPAR_OBJLIM:
      if ( lpi->spx->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
         *dval = lpi->spx->realParam(SoPlex::OBJLIMIT_UPPER);
      else
         *dval = lpi->spx->realParam(SoPlex::OBJLIMIT_LOWER);
      break;
   case SCIP_LPPAR_LPTILIM:
      *dval = lpi->spx->realParam(SoPlex::TIMELIMIT);
      break;
   case SCIP_LPPAR_ROWREPSWITCH:
      *dval = lpi->spx->realParam(SoPlex::REPRESENTATION_SWITCH);
      if( *dval >= SCIPlpiInfinity(lpi) )
         *dval = -1.0;
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      *dval = lpi->conditionlimit;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      lpi->spx->setFeastol(dval);
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      lpi->spx->setOpttol(dval);
      break;
   case SCIP_LPPAR_OBJLIM:
      if ( lpi->spx->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE )
         (void) lpi->spx->setRealParam(SoPlex::OBJLIMIT_UPPER, dval);
      else
         (void) lpi->spx->setRealParam(SoPlex::OBJLIMIT_LOWER, dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      (void) lpi->spx->setRealParam(SoPlex::TIMELIMIT, dval);
      break;
   case SCIP_LPPAR_ROWREPSWITCH:
      assert(dval >= -1.5);
      if( dval < 0.0 )
         (void) lpi->spx->setRealParam(SoPlex::REPRESENTATION_SWITCH, SCIPlpiInfinity(lpi));
      else
         (void) lpi->spx->setRealParam(SoPlex::REPRESENTATION_SWITCH, dval);
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      lpi->conditionlimit = dval;
      lpi->checkcondition = (dval >= 0);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);
   SCIPdebugMessage("calling SCIPlpiInfinity()\n");

   return lpi->spx->realParam(SoPlex::INFTY);
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val
   )
{
   assert(lpi != NULL);
   SCIPdebugMessage("calling SCIPlpiIsInfinity()\n");

   return (val >= lpi->spx->realParam(SoPlex::INFTY));
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** returns, whether the given file exists */
static
SCIP_Bool fileExists(
   const char*           filename            /**< file name */
   )
{
   FILE* f;

   f = fopen(filename, "r");
   if( f == NULL )
      return FALSE;

   fclose(f);

   return TRUE;
}

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   if( !fileExists(fname) )
      return SCIP_NOFILE;

   try
   {
      assert(lpi->spx->intParam(SoPlex::READMODE) == SoPlex::READMODE_REAL);
      if( !lpi->spx->readFile(fname) )
         return SCIP_READERROR;
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(fname != NULL);

   try
   {
      (void) lpi->spx->writeFileReal(fname);
   }
#ifndef NDEBUG
   catch( const SPxException& x )
   {
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#else
   catch( const SPxException& )
   {
#endif
      return SCIP_WRITEERROR;
   }

   return SCIP_OKAY;
}

/**@} */
