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

/**@file   reader_nl.cpp
 * @ingroup DEFPLUGINS_READER
 * @brief  AMPL .nl file reader
 * @author Stefan Vigerske
 *
 * For documentation on ampl::mp, see https://ampl.github.io and https://www.zverovich.net/2014/09/19/reading-nl-files.html.
 * For documentation on .nl files, see https://ampl.com/REFS/hooking2.pdf.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string>
#include <vector>
#include <map>

#include "scip/reader_nl.h"
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "scip/expr_sum.h"
#include "scip/expr_product.h"
#include "scip/expr_pow.h"
#include "scip/expr_log.h"
#include "scip/expr_exp.h"
#include "scip/expr_trig.h"
#include "scip/expr_abs.h"

// disable -Wshadow warnings for upcoming includes of AMPL/MP
// disable -Wimplicit-fallthrough as I don't want to maintain extra comments in AMPL/MP code to suppress these
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#if __GNUC__ >= 7
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
#endif

#include "mp/nl-reader.h"

#define READER_NAME             "nlreader"
#define READER_DESC             "AMPL .nl file reader"
#define READER_EXTENSION        "nl"

// a variant of SCIP_CALL that throws a std::logic_error if not SCIP_OKAY
// (using cast to long long to work around issues with old MSVC)
#define SCIP_CALL_THROW(x) \
   do                                                                                                   \
   {                                                                                                    \
      SCIP_RETCODE throw_retcode;                                                                       \
      if( ((throw_retcode) = (x)) != SCIP_OKAY )                                                        \
         throw std::logic_error("Error <" + std::to_string((long long)throw_retcode) + "> in function call"); \
   }                                                                                                    \
   while( false )

/*
 * Data structures
 */

/// problem data stored in SCIP
struct SCIP_ProbData
{
   char*                 filenamestub;       /**< name of input file, without .nl extension; array is long enough to hold 5 extra chars */
   int                   filenamestublen;    /**< length of filenamestub string */

   int                   amplopts[mp::MAX_AMPL_OPTIONS];  /**< AMPL options from .nl header */
   int                   namplopts;          /**< number of AMPL options from .nl header */

   SCIP_VAR**            vars;               /**< variables in the order given by AMPL */
   int                   nvars;              /**< number of variables */

   SCIP_CONS**           conss;              /**< constraints in the order given by AMPL */
   int                   nconss;             /**< number of constraints */

   SCIP_Bool             islp;               /**< whether problem is an LP (only linear constraints, only continuous vars) */
};

/*
 * Local methods
 */

// forward declaration
static SCIP_DECL_PROBDELORIG(probdataDelOrigNl);

/// implementation of AMPL/MPs NLHandler that constructs a SCIP problem while a .nl file is read
class AMPLProblemHandler : public mp::NLHandler<AMPLProblemHandler, SCIP_EXPR*>
{
private:
   SCIP* scip;
   SCIP_PROBDATA* probdata;

   // variable expressions corresponding to nonlinear variables
   // created in OnHeader() and released in destructor
   // for reuse of var-expressions in OnVariableRef()
   std::vector<SCIP_EXPR*> varexprs;

   // linear parts for nonlinear constraints
   // first collect and then add to constraints in EndInput()
   std::vector<std::vector<std::pair<SCIP_Real, SCIP_VAR*> > > nlconslin;

   // expression that represents a nonlinear objective function
   // used to create a corresponding constraint in EndInput(), unless NULL
   SCIP_EXPR* objexpr;

   // common expressions (defined variables from statements like "var xsqr = x^2;" in an AMPL model)
   // they are constructed by BeginCommonExpr/EndCommonExpr below and are referenced by index in OnCommonExprRef
   std::vector<SCIP_EXPR*> commonexprs;

   // collect expressions that need to be released eventually
   // this are all expression that are returned to the AMPL/MP code in AMPLProblemHandler::OnXyz() functions
   // they need to be released exactly once, but after they are used in another expression or a constraint
   // as AMPL/MP may reuse expressions (common subexpressions), we don't release an expression when it is used
   // as a child or when constructing a constraint, but first collect them all and then release in destructor
   // alternatively, one could encapsulate SCIP_EXPR* into a small class that handles proper reference counting
   std::vector<SCIP_EXPR*> exprstorelease;

   // SOS constraints
   // collected while handling suffixes in SuffixHandler
   // sosvars maps the SOS index (can be negative) to the indices of the variables in the SOS
   // sosweights gives for each variable its weight in the SOS it appears in (if any)
   std::map<int, std::vector<int> > sosvars;
   std::vector<int> sosweights;

   // initial solution, if any
   SCIP_SOL* initsol;

   // opened files with column/variable and row/constraint names, or NULL
   fmt::File* colfile;
   fmt::File* rowfile;

   // get name from names strings, if possible
   // returns whether a name has been stored
   bool nextName(
      const char*&       namesbegin,         /**< current pointer into names string, or NULL */
      const char*        namesend,           /**< pointer to end of names string */
      char*              name                /**< buffer to store name, should have length SCIP_MAXSTRLEN */
   )
   {
      if( namesbegin == NULL )
         return false;

      // copy namesbegin into name until newline or namesend
      // updates namesbegin
      int nchars = 0;
      while( namesbegin != namesend )
      {
         if( nchars == SCIP_MAXSTRLEN )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "name too long when parsing names file");
            // do no longer read names from this string (something seems awkward)
            namesbegin = NULL;
            return false;
         }
         if( *namesbegin == '\n' )
         {
            *name = '\0';
            ++namesbegin;
            return true;
         }
         *(name++) = *(namesbegin++);
         ++nchars;
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "missing newline when parsing names file");
      return false;
   }

public:
   /// constructor
   ///
   /// initializes SCIP problem and problem data
   AMPLProblemHandler(
      SCIP*              scip_,              ///< SCIP data structure
      const char*        filename            ///< name of .nl file that is read
      )
   : scip(scip_),
     probdata(NULL),
     objexpr(NULL),
     initsol(NULL),
     colfile(NULL),
     rowfile(NULL)
   {
      assert(scip != NULL);
      assert(filename != NULL);

      SCIP_CALL_THROW( SCIPallocClearMemory(scip, &probdata) );

      /* get name of input file without file extension (if any) */
      const char* extstart = strrchr(const_cast<char*>(filename), '.');
      if( extstart != NULL )
         probdata->filenamestublen = extstart - filename;
      else
         probdata->filenamestublen = strlen(filename);
      assert(probdata->filenamestublen > 0);
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &probdata->filenamestub, probdata->filenamestublen + 5) );
      memcpy(probdata->filenamestub, filename, probdata->filenamestublen);
      probdata->filenamestub[probdata->filenamestublen] = '\0';

      /* derive probname from name of input file without path and extension */
      const char* probname = strrchr(probdata->filenamestub, '/');
      if( probname == NULL )
         probname = probdata->filenamestub;
      else
         ++probname;

      // initialize empty SCIP problem
      SCIP_CALL_THROW( SCIPcreateProb(scip, probname, probdataDelOrigNl, NULL, NULL, NULL, NULL, NULL, probdata) );

      // try to open files with variable and constraint names
      // temporarily add ".col" and ".row", respectively, to filenamestub
      try
      {
         probdata->filenamestub[probdata->filenamestublen] = '.';
         probdata->filenamestub[probdata->filenamestublen+1] = 'c';
         probdata->filenamestub[probdata->filenamestublen+2] = 'o';
         probdata->filenamestub[probdata->filenamestublen+3] = 'l';
         probdata->filenamestub[probdata->filenamestublen+4] = '\0';
         colfile = new fmt::File(probdata->filenamestub, fmt::File::RDONLY);

         probdata->filenamestub[probdata->filenamestublen+1] = 'r';
         probdata->filenamestub[probdata->filenamestublen+3] = 'w';
         rowfile = new fmt::File(probdata->filenamestub, fmt::File::RDONLY);
      }
      catch( const fmt::SystemError& e )
      {
         // probably a file open error, probably because file not found
         // ignore, we can make up our own names
      }
      probdata->filenamestub[probdata->filenamestublen] = '\0';
   }

   AMPLProblemHandler(const AMPLProblemHandler&) = delete;
   AMPLProblemHandler& operator=(const AMPLProblemHandler&) = delete;

   /// destructor
   ///
   /// only asserts that cleanup() has been called, as we cannot throw an exception or return a SCIP_RETCODE here
   ~AMPLProblemHandler()
   {
      // exprs and linear constraint arrays should have been cleared up in cleanup()
      assert(varexprs.empty());
      assert(exprstorelease.empty());

      delete colfile;
      delete rowfile;
   }

   /// process header of .nl files
   ///
   /// create and add variables, allocate constraints
   void OnHeader(
      const mp::NLHeader& h                  ///< header data
      )
   {
      char name[SCIP_MAXSTRLEN];
      int nnlvars;

      assert(probdata->vars == NULL);
      assert(probdata->conss == NULL);

      probdata->namplopts = h.num_ampl_options;
      BMScopyMemoryArray(probdata->amplopts, h.ampl_options, h.num_ampl_options);

      // read variable and constraint names from file, if available, into memory
      // if not available, we will get varnamesbegin==NULL and consnamesbegin==NULL
      mp::MemoryMappedFile<> mapped_colfile;
      if( colfile != NULL )
         mapped_colfile.map(*colfile, "colfile");
      const char* varnamesbegin = mapped_colfile.start();
      const char* varnamesend = mapped_colfile.start() + mapped_colfile.size();

      mp::MemoryMappedFile<> mapped_rowfile;
      if( rowfile != NULL )
         mapped_rowfile.map(*rowfile, "rowfile");
      const char* consnamesbegin = mapped_rowfile.start();
      const char* consnamesend = mapped_rowfile.start() + mapped_rowfile.size();

      probdata->nvars = h.num_vars;
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &probdata->vars, probdata->nvars) );

      // number of nonlinear variables
      nnlvars = MAX(h.num_nl_vars_in_cons, h.num_nl_vars_in_objs);
      varexprs.resize(nnlvars);

      // create variables
      // create variable expressions for nonlinear variables
      for( int i = 0; i < h.num_vars; ++i )
      {
         SCIP_VARTYPE vartype;
         // Nonlinear variables in both constraints and objective
         if( i < h.num_nl_vars_in_both - h.num_nl_integer_vars_in_both )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_nl_vars_in_both )
            vartype = SCIP_VARTYPE_INTEGER;
         // Nonlinear variables in constraints
         else if( i < h.num_nl_vars_in_cons - h.num_nl_integer_vars_in_cons )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_nl_vars_in_cons )
            vartype = SCIP_VARTYPE_INTEGER;
         // Nonlinear variables in objective
         else if( i < h.num_nl_vars_in_objs - h.num_nl_integer_vars_in_objs )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_nl_vars_in_objs )
            vartype = SCIP_VARTYPE_INTEGER;
         // Linear variables
         else if( i < h.num_vars - h.num_linear_binary_vars - h.num_linear_integer_vars )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_vars - h.num_linear_integer_vars )
            vartype = SCIP_VARTYPE_BINARY;
         else
            vartype = SCIP_VARTYPE_INTEGER;

         if( !nextName(varnamesbegin, varnamesend, name) )
         {
            // make up name if no names file or could not be read
            switch( vartype )
            {
               case SCIP_VARTYPE_BINARY :
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b%d", i);
                  break;
               case SCIP_VARTYPE_INTEGER :
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "i%d", i);
                  break;
               case SCIP_VARTYPE_CONTINUOUS :
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
                  break;
               // coverity[deadcode]
               default:
                  SCIPABORT();
                  break;
            }
         }

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &probdata->vars[i], name,
            vartype == SCIP_VARTYPE_BINARY ? 0.0 : -SCIPinfinity(scip),
            vartype == SCIP_VARTYPE_BINARY ? 1.0 :  SCIPinfinity(scip),
            0.0, vartype) );
         SCIP_CALL_THROW( SCIPaddVar(scip, probdata->vars[i]) );

         if( i < nnlvars )
         {
            SCIP_CALL_THROW( SCIPcreateExprVar(scip, &varexprs[i], probdata->vars[i], NULL, NULL) );
         }
      }

      // alloc some space for constraints
      probdata->nconss = h.num_algebraic_cons;
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &probdata->conss, probdata->nconss) );
      nlconslin.resize(h.num_nl_cons);

      // create empty nonlinear constraints
      // use expression == 0, because nonlinear constraint don't like to be without an expression
      SCIP_EXPR* dummyexpr;
      SCIP_CALL_THROW( SCIPcreateExprValue(scip, &dummyexpr, 0.0, NULL, NULL) );
      for( int i = 0; i < h.num_nl_cons; ++i )
      {
         // make up name if no names file or could not be read
         if( !nextName(consnamesbegin, consnamesend, name) )
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlc%d", i);

         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &probdata->conss[i], name, dummyexpr, -SCIPinfinity(scip), SCIPinfinity(scip)) );
      }
      SCIP_CALL_THROW( SCIPreleaseExpr(scip, &dummyexpr) );

      // create empty linear constraints
      for( int i = h.num_nl_cons; i < h.num_algebraic_cons; ++i )
      {
         if( !nextName(consnamesbegin, consnamesend, name) )
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lc%d", i);
         SCIP_CALL_THROW( SCIPcreateConsBasicLinear(scip, &probdata->conss[i], name, 0, NULL, NULL, -SCIPinfinity(scip), SCIPinfinity(scip)) );
      }

      if( h.num_nl_cons == 0 && h.num_integer_vars() == 0 )
         probdata->islp = true;

      // alloc space for common expressions
      commonexprs.resize(h.num_common_exprs());
   }

   /// receive notification of a number in a nonlinear expression
   SCIP_EXPR* OnNumber(
      double             value               ///< value
      )
   {
      SCIP_EXPR* expr;

      SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, value, NULL, NULL) );

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// receive notification of a variable reference in a nonlinear expression
   SCIP_EXPR* OnVariableRef(
      int                variableIndex       ///< AMPL index of variable
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < (int)varexprs.size());
      assert(varexprs[variableIndex] != NULL);

      return varexprs[variableIndex];
   }

   /// receive notification of a unary expression
   SCIP_EXPR* OnUnary(
      mp::expr::Kind     kind,               ///< expression operator
      SCIP_EXPR*         child               ///< argument
      )
   {
      SCIP_EXPR* expr;

      assert(child != NULL);

      switch( kind )
      {
         case mp::expr::MINUS:
         {
            SCIP_Real minusone = -1.0;
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 1, &child, &minusone, 0.0, NULL, NULL) );
            break;
         }

         case mp::expr::ABS:
            SCIP_CALL_THROW( SCIPcreateExprAbs(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::POW2:
            SCIP_CALL_THROW( SCIPcreateExprPow(scip, &expr, child, 2.0, NULL, NULL) );
            break;

         case mp::expr::SQRT:
            SCIP_CALL_THROW( SCIPcreateExprPow(scip, &expr, child, 0.5, NULL, NULL) );
            break;

         case mp::expr::LOG:
            SCIP_CALL_THROW( SCIPcreateExprLog(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::LOG10:  // 1/log(10)*log(child)
         {
            SCIP_EXPR* logexpr;
            SCIP_Real factor = 1.0/log(10.0);
            SCIP_CALL_THROW( SCIPcreateExprLog(scip, &logexpr, child, NULL, NULL) );
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 1, &logexpr, &factor, 0.0, NULL, NULL) );
            SCIP_CALL_THROW( SCIPreleaseExpr(scip, &logexpr) );
            break;
         }

         case mp::expr::EXP:
            SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::SIN:
            SCIP_CALL_THROW( SCIPcreateExprSin(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::COS:
            SCIP_CALL_THROW( SCIPcreateExprCos(scip, &expr, child, NULL, NULL) );
            break;

         default:
            OnUnhandled(mp::expr::str(kind));
            return NULL;
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// receive notification of a binary expression
   SCIP_EXPR* OnBinary(
      mp::expr::Kind     kind,               ///< expression operand
      SCIP_EXPR*         firstChild,         ///< first argument
      SCIP_EXPR*         secondChild         ///< second argument
      )
   {
      SCIP_EXPR* expr;
      SCIP_EXPR* children[2] = { firstChild, secondChild };

      assert(firstChild != NULL);
      assert(secondChild != NULL);

      switch( kind )
      {
         case mp::expr::ADD:
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 2, children, NULL, 0.0, NULL, NULL) );
            break;

         case mp::expr::SUB:
         {
            SCIP_Real coefs[2] = { 1.0, -1.0 };
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 2, children, coefs, 0.0, NULL, NULL) );
            break;
         }

         case mp::expr::MUL:
            SCIP_CALL_THROW( SCIPcreateExprProduct(scip, &expr, 2, children, 1.0, NULL, NULL) );
            break;

         case mp::expr::DIV:
            SCIP_CALL_THROW( SCIPcreateExprPow(scip, &children[1], secondChild, -1.0, NULL, NULL) );
            SCIP_CALL_THROW( SCIPcreateExprProduct(scip, &expr, 2, children, 1.0, NULL, NULL) );
            SCIP_CALL_THROW( SCIPreleaseExpr(scip, &children[1]) );
            break;

         case mp::expr::POW_CONST_BASE:
         case mp::expr::POW_CONST_EXP:
         case mp::expr::POW:
            // with some .nl files, we seem to get mp::expr::POW even if base or exponent is constant,
            // so do not rely on kind but better check expr type
            if( SCIPisExprValue(scip, secondChild) )
            {
               SCIP_CALL_THROW( SCIPcreateExprPow(scip, &expr, firstChild, SCIPgetValueExprValue(secondChild), NULL, NULL) );
               break;
            }

            if( SCIPisExprValue(scip, firstChild) && SCIPgetValueExprValue(firstChild) > 0.0 )
            {
               // reformulate constant^y as exp(y*log(constant)), if constant > 0.0
               // if constant < 0, we create an expression and let cons_nonlinear figure out infeasibility somehow
               SCIP_EXPR* prod;

               SCIP_Real coef = log(SCIPgetValueExprValue(firstChild)); // log(firstChild)
               SCIP_CALL_THROW( SCIPcreateExprSum(scip, &prod, 1, &secondChild, &coef, 0.0, NULL, NULL) );  // log(firstChild)*secondChild
               SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, prod, NULL, NULL) );  // expr(log(firstChild)*secondChild)

               SCIP_CALL_THROW( SCIPreleaseExpr(scip, &prod) );
               break;
            }

            {
               // reformulate x^y as exp(y*log(x))
               SCIP_EXPR* prod;

               assert(SCIPisExprValue(scip, secondChild));

               SCIP_CALL_THROW( SCIPcreateExprLog(scip, &children[0], firstChild, NULL, NULL) );  // log(firstChild)
               SCIP_CALL_THROW( SCIPcreateExprProduct(scip, &prod, 2, children, 1.0, NULL, NULL) );  // log(firstChild)*secondChild
               SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, prod, NULL, NULL) );  // expr(log(firstChild)*secondChild)

               SCIP_CALL_THROW( SCIPreleaseExpr(scip, &prod) );
               SCIP_CALL_THROW( SCIPreleaseExpr(scip, &children[0]) );
               break;
            }

         default:
            OnUnhandled(mp::expr::str(kind));
            return NULL;
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// handler to create a list of terms in a sum
   ///
   /// NumericArgHandler is copied around, so it keeps only a pointer (with reference counting) to actual data
   class NumericArgHandler
   {
   public:
      std::shared_ptr<std::vector<SCIP_EXPR*> > v;

      /// constructor
      NumericArgHandler(
         int             num_args            ///< number of terms to expect
         )
      : v(new std::vector<SCIP_EXPR*>())
      {
         v->reserve(num_args);
      }

      /// adds term to sum
      void AddArg(
         SCIP_EXPR*      term                ///< term to add
         )
      {
         v->push_back(term);
      }
   };

   /// receive notification of the beginning of a summation
   NumericArgHandler BeginSum(
      int                num_args            ///< number of terms to expect
      )
   {
      NumericArgHandler h(num_args);
      return h;
   }

   /// receive notification of the end of a summation
   SCIP_EXPR* EndSum(
      NumericArgHandler  handler             ///< handler that handled the sum
      )
   {
      SCIP_EXPR* expr;
      SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, (int)handler.v->size(), handler.v->data(), NULL, 0.0, NULL, NULL) );
      // remember that we have to release this expr
      exprstorelease.push_back(expr);
      return expr;
   }

   /// receive notification of an objective type and the nonlinear part of an objective expression
   void OnObj(
      int                objectiveIndex,     ///< index of objective
      mp::obj::Type      type,               ///< objective sense
      SCIP_EXPR*         nonlinearExpression ///< nonlinear part of objective function
      )
   {
      if( objectiveIndex >= 1 )
         OnUnhandled("multiple objective functions");

      SCIP_CALL_THROW( SCIPsetObjsense(scip, type == mp::obj::Type::MAX ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE) );

      assert(objexpr == NULL);

      if( nonlinearExpression != NULL && SCIPisExprValue(scip, nonlinearExpression) )
      {
         // handle objective constant by adding a fixed variable for it
         SCIP_VAR* objconstvar;
         SCIP_Real objconst = SCIPgetValueExprValue(nonlinearExpression);

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &objconstvar, "objconstant", objconst, objconst, 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL_THROW( SCIPaddVar(scip, objconstvar) );
         SCIP_CALL_THROW( SCIPreleaseVar(scip, &objconstvar) );
      }
      else
      {
         objexpr = nonlinearExpression;
      }
   }

   /// receive notification of an algebraic constraint expression
   void OnAlgebraicCon(
      int                constraintIndex,    ///< index of constraint
      SCIP_EXPR*         expr                ///< nonlinear part of constraint
      )
   {
      if( expr != NULL )
      {
         SCIP_CALL_THROW( SCIPchgExprNonlinear(scip, probdata->conss[constraintIndex], expr) );
      }
   }

   /// handles linear part of a common expression
   /// sets up a sum expression, if the linear part isn't empty
   class LinearExprHandler
   {
   private:
      AMPLProblemHandler& amplph;
      SCIP_EXPR*          commonexpr;

   public:
      /// constructor
      LinearExprHandler(
         AMPLProblemHandler& amplph_,        ///< problem handler
         int                 index,          ///< index of common expression
         int                 num_linear_terms///< number of terms to expect
         )
      : amplph(amplph_),
        commonexpr(NULL)
      {
         if( num_linear_terms > 0 )
         {
            SCIP_CALL_THROW( SCIPcreateExprSum(amplph.scip, &commonexpr, 0, NULL, NULL, 0.0, NULL, NULL) );
            amplph.commonexprs[index] = commonexpr;
            amplph.exprstorelease.push_back(commonexpr);
         }
      }

      /// receives notification of a term in the linear expression
      void AddTerm(
         int             var_index,          ///< AMPL index of variable
         double          coef                ///< variable coefficient
         )
      {
         assert(commonexpr != NULL);

         if( coef == 0.0 )
            return;

         if( var_index < (int)amplph.varexprs.size() )
         {
            SCIP_CALL_THROW( SCIPappendExprSumExpr(amplph.scip, commonexpr, amplph.varexprs[var_index], coef) );
         }
         else
         {
            // the index variable is linear (not sure this can happen here)
            assert(var_index < amplph.probdata->nvars);
            SCIP_EXPR* varexpr;
            SCIP_CALL_THROW( SCIPcreateExprVar(amplph.scip, &varexpr, amplph.probdata->vars[var_index], NULL, NULL) );
            SCIP_CALL_THROW( SCIPappendExprSumExpr(amplph.scip, commonexpr, varexpr, coef) );
            SCIP_CALL_THROW( SCIPreleaseExpr(amplph.scip, &varexpr) );
         }
      }
   };

   /// receive notification of the beginning of a common expression (defined variable)
   LinearExprHandler BeginCommonExpr(
      int                index,              ///< index of common expression
      int                num_linear_terms    ///< number of terms to expect
      )
   {
      assert(index >= 0);
      assert(index < (int)commonexprs.size());

      return LinearExprHandler(*this, index, num_linear_terms);
   }

   /// receive notification of the end of a common expression
   void EndCommonExpr(
      int                index,              ///< index of common expression
      SCIP_EXPR*         expr,               ///< nonlinear part of common expression
      int                /* position */      ///< argument that doesn't seem to have any purpose
      )
   {
      if( commonexprs[index] != NULL )
      {
         // add expr, if any, to linear part
         if( expr != NULL )
         {
            SCIP_CALL_THROW( SCIPappendExprSumExpr(scip, commonexprs[index], expr, 1.0) );
         }
      }
      else if( expr != NULL )
      {
         commonexprs[index] = expr;
      }
   }

   /// receive notification of a common expression (defined variable) reference
   SCIP_EXPR* OnCommonExprRef(
      int                expr_index          ///< index of common expression
      )
   {
      assert(expr_index >= 0);
      assert(expr_index < (int)commonexprs.size());
      assert(commonexprs[expr_index] != NULL);
      return commonexprs[expr_index];
   }

   /// receive notification of variable bounds
   void OnVarBounds(
      int                variableIndex,      ///< AMPL index of variable
      double             variableLB,         ///< variable lower bound
      double             variableUB          ///< variable upper bound
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < probdata->nvars);

      // as far as I see, ampl::mp gives -inf, +inf for no-bounds, which is always beyond SCIPinfinity()
      // we ignore bounds outside [-scipinfinity,scipinfinity] here
      // for binary variables, we also ignore bounds outside [0,1]
      if( variableLB > (SCIPvarGetType(probdata->vars[variableIndex]) == SCIP_VARTYPE_BINARY ? 0.0 : -SCIPinfinity(scip)) )
      {
         SCIP_CALL_THROW( SCIPchgVarLbGlobal(scip, probdata->vars[variableIndex], variableLB) );
      }
      if( variableUB < (SCIPvarGetType(probdata->vars[variableIndex]) == SCIP_VARTYPE_BINARY ? 1.0 :  SCIPinfinity(scip)) )
      {
         SCIP_CALL_THROW( SCIPchgVarUbGlobal(scip, probdata->vars[variableIndex], variableUB) );
      }
   }

   /// receive notification of constraint sides
   void OnConBounds(
      int                index,              ///< AMPL index of constraint
      double             lb,                 ///< constraint left-hand-side
      double             ub                  ///< constraint right-hand-side
      )
   {
      assert(index >= 0);
      assert(index < probdata->nconss);

      // nonlinear constraints are first
      if( index < (int)nlconslin.size() )
      {
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_CALL_THROW( SCIPchgLhsNonlinear(scip, probdata->conss[index], lb) );
         }
         if( !SCIPisInfinity(scip,  ub) )
         {
            SCIP_CALL_THROW( SCIPchgRhsNonlinear(scip, probdata->conss[index], ub) );
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_CALL_THROW( SCIPchgLhsLinear(scip, probdata->conss[index], lb) );
         }
         if( !SCIPisInfinity(scip,  ub) )
         {
            SCIP_CALL_THROW( SCIPchgRhsLinear(scip, probdata->conss[index], ub) );
         }
      }
   }

   /// receive notification of the initial value for a variable
   void OnInitialValue(
      int                var_index,          ///< AMPL index of variable
      double             value               ///< initial primal value of variable
      )
   {
      if( initsol == NULL )
      {
         SCIP_CALL_THROW( SCIPcreateSol(scip, &initsol, NULL) );
      }

      SCIP_CALL_THROW( SCIPsetSolVal(scip, initsol, probdata->vars[var_index], value) );
   }

   /// receives notification of the initial value for a dual variable
   void OnInitialDualValue(
      int                /* con_index */,    ///< AMPL index of constraint
      double             /* value */         ///< initial dual value of constraint
      )
   {
      // ignore initial dual value
   }

   /// receives notification of Jacobian column sizes
   ColumnSizeHandler OnColumnSizes()
   {
      /// use ColumnSizeHandler from upper class, which does nothing
      return ColumnSizeHandler();
   }

   /// handling of suffices for variable and constraint flags and SOS constraints
   ///
   /// regarding SOS in AMPL, see https://ampl.com/faqs/how-can-i-use-the-solvers-special-ordered-sets-feature/
   /// we pass the .ref suffix as weight to the SOS constraint handlers
   /// for a SOS2, the weights determine the order of variables in the set
   template<typename T> class SuffixHandler
   {
   private:
      AMPLProblemHandler& amplph;

      // type of suffix that is handled, or IGNORE if unsupported suffix
      enum
      {
         IGNORE,
         CONSINITIAL,
         CONSSEPARATE,
         CONSENFORCE,
         CONSCHECK,
         CONSPROPAGATE,
         CONSDYNAMIC,
         CONSREMOVABLE,
         VARINITIAL,
         VARREMOVABLE,
         VARSOSNO,
         VARREF,
      } suffix;

   public:
      /// constructor
      SuffixHandler(
         AMPLProblemHandler& amplph_,        ///< problem handler
         fmt::StringRef      name,           ///< name of suffix
         mp::suf::Kind       kind            ///< whether suffix applies to var, cons, etc
         )
      : amplph(amplph_),
        suffix(IGNORE)
      {
         switch( kind )
         {
            case mp::suf::Kind::CON:
               if( strncmp(name.data(), "initial", name.size()) == 0 )
               {
                  suffix = CONSINITIAL;
               }
               else if( strncmp(name.data(), "separate", name.size()) == 0 )
               {
                  suffix = CONSSEPARATE;
               }
               else if( strncmp(name.data(), "enforce", name.size()) == 0 )
               {
                  suffix = CONSENFORCE;
               }
               else if( strncmp(name.data(), "check", name.size()) == 0 )
               {
                  suffix = CONSCHECK;
               }
               else if( strncmp(name.data(), "propagate", name.size()) == 0 )
               {
                  suffix = CONSPROPAGATE;
               }
               else if( strncmp(name.data(), "dynamic", name.size()) == 0 )
               {
                  suffix = CONSDYNAMIC;
               }
               else if( strncmp(name.data(), "removable", name.size()) == 0 )
               {
                  suffix = CONSREMOVABLE;
               }
               else
               {
                  SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown constraint suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               }
               break;

            case mp::suf::Kind::VAR:
            {
               if( strncmp(name.data(), "initial", name.size()) == 0 )
               {
                  suffix = VARINITIAL;
               }
               else if( strncmp(name.data(), "removable", name.size()) == 0 )
               {
                  suffix = VARREMOVABLE;
               }
               else if( strncmp(name.data(), "sosno", name.size()) == 0 )
               {
                  // SOS membership
                  suffix = VARSOSNO;
               }
               else if( strncmp(name.data(), "ref", name.size()) == 0 )
               {
                  // SOS weights
                  suffix = VARREF;
                  amplph.sosweights.resize(amplph.probdata->nvars, 0);
               }
               else
               {
                  SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown variable suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               }
               break;

            case mp::suf::Kind::OBJ:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown objective suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;

            case mp::suf::Kind::PROBLEM:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown problem suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;
            }
         }
      }

      void SetValue(
         int             index,              ///< index of variable, constraint, etc
         T               value               ///< value of suffix
      )
      {
         assert(index >= 0);
         switch( suffix )
         {
            case IGNORE :
               return;

            case CONSINITIAL:
               SCIP_CALL_THROW( SCIPsetConsInitial(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSSEPARATE:
               SCIP_CALL_THROW( SCIPsetConsSeparated(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSENFORCE:
               SCIP_CALL_THROW( SCIPsetConsEnforced(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSCHECK:
               SCIP_CALL_THROW( SCIPsetConsChecked(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSPROPAGATE:
               SCIP_CALL_THROW( SCIPsetConsPropagated(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSDYNAMIC:
               SCIP_CALL_THROW( SCIPsetConsDynamic(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSREMOVABLE:
               SCIP_CALL_THROW( SCIPsetConsRemovable(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case VARINITIAL:
               assert(index < amplph.probdata->nvars);
               SCIP_CALL_THROW( SCIPvarSetInitial(amplph.probdata->vars[index], value == 1) );
               break;

            case VARREMOVABLE:
               assert(index < amplph.probdata->nvars);
               SCIP_CALL_THROW( SCIPvarSetRemovable(amplph.probdata->vars[index], value == 1) );
               break;

            case VARSOSNO:
               // remember that variable index belongs to SOS identified by value
               amplph.sosvars[(int)value].push_back(index);
               break;

            case VARREF:
               // remember that variable index has weight value
               amplph.sosweights[index] = (int)value;
               break;
         }
      }
   };

   typedef SuffixHandler<int> IntSuffixHandler;
   /// receive notification of an integer suffix
   IntSuffixHandler OnIntSuffix(
      fmt::StringRef     name,               ///< suffix name, not null-terminated
      mp::suf::Kind      kind,               ///< suffix kind
      int                /*num_values*/      ///< number of values to expect
      )
   {
      return IntSuffixHandler(*this, name, kind);
   }

   typedef SuffixHandler<SCIP_Real> DblSuffixHandler;
   /// receive notification of a double suffix
   DblSuffixHandler OnDblSuffix(
      fmt::StringRef     name,               ///< suffix name, not null-terminated
      mp::suf::Kind      kind,               ///< suffix kind
      int                /*num_values*/      ///< number of values to expect
      )
   {
      return DblSuffixHandler(*this, name, kind);
   }

   /// handles receiving the linear part of an objective or constraint
   ///
   /// for objective, set the objective-coefficient of the variable
   /// for linear constraints, add to the constraint
   /// for nonlinear constraints, add to nlconslin vector; adding to constraint later
   class LinearPartHandler
   {
   private:
      AMPLProblemHandler& amplph;
      int constraintIndex;

   public:
      // constructor for constraint
      explicit LinearPartHandler(
         AMPLProblemHandler& amplph_,        ///< problem handler
         int                 constraintIndex_///< constraint index
         )
      : amplph(amplph_),
        constraintIndex(constraintIndex_)
      {
         assert(constraintIndex_ >= 0);
         assert(constraintIndex_ < amplph.probdata->nconss);
      }

      // constructor for linear objective
      explicit LinearPartHandler(
         AMPLProblemHandler& amplph_         ///< problem handler
         )
      : amplph(amplph_),
        constraintIndex(-1)
      { }

      void AddTerm(
         int             variableIndex,      ///< AMPL index of variable
         double          coefficient         ///< coefficient of variable
         )
      {
         assert(variableIndex >= 0);
         assert(variableIndex < amplph.probdata->nvars);

         if( coefficient == 0.0 )
            return;

         if( constraintIndex < 0 )
         {
            SCIP_CALL_THROW( SCIPchgVarObj(amplph.scip, amplph.probdata->vars[variableIndex], coefficient) );
         }
         else if( constraintIndex < (int)amplph.nlconslin.size() )
         {
            amplph.nlconslin[constraintIndex].push_back(std::pair<SCIP_Real, SCIP_VAR*>(coefficient, amplph.probdata->vars[variableIndex]));
         }
         else
         {
            SCIP_CONS* lincons = amplph.probdata->conss[constraintIndex];
            SCIP_CALL_THROW( SCIPaddCoefLinear(amplph.scip, lincons, amplph.probdata->vars[variableIndex], coefficient) );
         }
      }
   };

   typedef LinearPartHandler LinearObjHandler;

   /// receive notification of the linear part of an objective
   LinearPartHandler OnLinearObjExpr(
      int                objectiveIndex,     ///< index of objective
      int                /* numLinearTerms *////< number of terms to expect
      )
   {
      if( objectiveIndex >= 1 )
         OnUnhandled("multiple objective functions");

      return LinearObjHandler(*this);
   }

   typedef LinearPartHandler LinearConHandler;

   /// receive notification of the linear part of a constraint
   LinearConHandler OnLinearConExpr(
      int                constraintIndex,    ///< index of constraint
      int                /* numLinearTerms *////< number of terms to expect
      )
   {
      return LinearConHandler(*this, constraintIndex);
   }

   /// receive notification of the end of the input
   ///
   /// - setup all nonlinear constraints and add them to SCIP
   /// - add linear constraints to SCIP (should be after nonlinear ones to respect order in .nl file)
   /// - add initial solution, if initial values were given
   void EndInput()
   {
      // turn nonlinear objective into constraint
      // min f(x) -> min z s.t. f(x) - z <= 0
      // max f(x) -> max z s.t. 0 <= f(x) - z
      if( objexpr != NULL )
      {
         SCIP_CONS* objcons;
         SCIP_VAR* objvar;

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &objvar, "nlobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL_THROW( SCIPaddVar(scip, objvar) );

         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &objcons, "objcons", objexpr,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0) );
         SCIP_CALL_THROW( SCIPaddLinearVarNonlinear(scip, objcons, objvar, -1.0) );
         SCIP_CALL_THROW( SCIPaddCons(scip, objcons) );

         SCIP_CALL_THROW( SCIPreleaseCons(scip, &objcons) );
         SCIP_CALL_THROW( SCIPreleaseVar(scip, &objvar) );
      }

      // add linear terms to expressions of nonlinear constraints (should be ok to do this one-by-one for now)
      for( size_t i = 0; i < nlconslin.size(); ++i )
      {
         for( size_t j = 0; j < nlconslin[i].size(); ++j )
         {
            SCIP_CALL_THROW( SCIPaddLinearVarNonlinear(scip, probdata->conss[i], nlconslin[i][j].second, nlconslin[i][j].first) );
         }
      }

      // add constraints
      for( int i = 0; i < probdata->nconss; ++i )
      {
         SCIP_CALL_THROW( SCIPaddCons(scip, probdata->conss[i]) );
      }

      // add SOS constraints
      std::vector<SCIP_VAR*> setvars;     // variables in one SOS
      std::vector<SCIP_Real> setweights;  // weights for one SOS
      if( !sosvars.empty() )
      {
         setvars.resize(probdata->nvars);
         probdata->islp = false;
      }
      if( !sosweights.empty() )
         setweights.resize(probdata->nvars);
      for( std::map<int, std::vector<int> >::iterator sosit(sosvars.begin()); sosit != sosvars.end(); ++sosit )
      {
         assert(sosit->first != 0);
         assert(!sosit->second.empty());

         // a negative SOS identifier means SOS2
         bool issos2 = sosit->first < 0;

         if( issos2 && sosweights.empty() )
         {
            // if no .ref suffix was given for a SOS2 constraint, then we consider this as an error
            // since the weights determine the order
            // for a SOS1, the weights only specify branching preference, so can treat them as optional
            OnUnhandled("SOS2 requires variable .ref suffix");
         }

         for( size_t i = 0; i < sosit->second.size(); ++i )
         {
            int varidx = sosit->second[i];
            setvars[i] = probdata->vars[varidx];

            if( issos2 && sosweights[varidx] == 0 )
               // 0 is the default if no ref was given for a variable; we don't allow this for SOS2
               OnUnhandled("Missing .ref value for SOS2 variable");
            if( !sosweights.empty() )
               setweights[i] = (SCIP_Real)sosweights[varidx];
         }

         SCIP_CONS* cons;
         char name[20];
         if( !issos2 )
         {
            (void) SCIPsnprintf(name, 20, "sos1_%d", sosit->first);
            SCIP_CALL_THROW( SCIPcreateConsBasicSOS1(scip, &cons, name, sosit->second.size(), setvars.data(), setweights.empty() ? NULL : setweights.data()) );
         }
         else
         {
            (void) SCIPsnprintf(name, 20, "sos2_%d", -sosit->first);
            SCIP_CALL_THROW( SCIPcreateConsBasicSOS2(scip, &cons, name, sosit->second.size(), setvars.data(), setweights.data()) );
         }
         SCIP_CALL_THROW( SCIPaddCons(scip, cons) );
         SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );
      }

      // add initial solution
      if( initsol != NULL )
      {
         SCIP_Bool stored;
         SCIP_CALL_THROW( SCIPaddSolFree(scip, &initsol, &stored) );
      }

      // release expressions
      SCIP_CALL_THROW( cleanup() );
   }

   /// releases expressions and linear constraints from data
   ///
   /// should be called if there was an error while reading the .nl file
   /// this is not in the destructor, because we want to return SCIP_RETCODE
   SCIP_RETCODE cleanup()
   {
      // release initial sol (in case EndInput() wasn't called)
      if( initsol != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &initsol) );
      }

      // release created expressions (they should all be used in other expressions or constraints now)
      while( !exprstorelease.empty() )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &exprstorelease.back()) );
         exprstorelease.pop_back();
      }

      // release variable expressions (they should all be used in other expressions or constraints now)
      while( !varexprs.empty() )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &varexprs.back()) );
         varexprs.pop_back();
      }

      return SCIP_OKAY;
   }
};


/*
 * Callback methods of probdata
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdataDelOrigNl)
{
   int i;

   assert((*probdata)->vars != NULL || (*probdata)->nvars == 0);
   assert((*probdata)->conss != NULL || (*probdata)->conss == 0);

   for( i = 0; i < (*probdata)->nconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->conss, (*probdata)->nconss);

   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->vars, (*probdata)->nvars);

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->filenamestub, (*probdata)->filenamestublen+5);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyNl)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderNl(scip) );

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadNl)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   try
   {
      // try to read the .nl file and setup SCIP problem
      AMPLProblemHandler handler(scip, filename);
      try
      {
         mp::ReadNLFile(filename, handler);
      }
      catch( const mp::UnsupportedError& e )
      {
         SCIPerrorMessage("unsupported construct in AMPL .nl file %s: %s\n", filename, e.what());

         SCIP_CALL( handler.cleanup() );

         return SCIP_READERROR;
      }
      catch( const mp::Error& e )
      {
         // some other error from ampl/mp, maybe invalid .nl file
         SCIPerrorMessage("%s\n", e.what());

         SCIP_CALL( handler.cleanup() );

         return SCIP_READERROR;
      }
      catch( const fmt::SystemError& e )
      {
         // probably a file open error, probably because file not found
         SCIPerrorMessage("%s\n", e.what());

         SCIP_CALL( handler.cleanup() );

         return SCIP_NOFILE;
      }
      catch( const std::bad_alloc& e )
      {
         SCIPerrorMessage("Out of memory: %s\n", e.what());

         SCIP_CALL( handler.cleanup() );

         return SCIP_NOMEMORY;
      }
   }
   catch( const std::exception& e )
   {
      SCIPerrorMessage("%s\n", e.what());
      return SCIP_ERROR;
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the AMPL .nl file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderNl(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_READER* reader = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );
   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyNl) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadNl) );

   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, "AMPL/MP 4e2d45c4", "AMPL .nl file reader library (github.com/ampl/mp)") );

   return SCIP_OKAY;
}

/** writes AMPL solution file
 *
 * problem must have been read with .nl reader
 */
SCIP_RETCODE SCIPwriteSolutionNl(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   if( probdata == NULL )
   {
      SCIPerrorMessage("No AMPL nl file read. Cannot write AMPL solution.\n");
      return SCIP_ERROR;
   }

   probdata->filenamestub[probdata->filenamestublen] = '.';
   probdata->filenamestub[probdata->filenamestublen+1] = 's';
   probdata->filenamestub[probdata->filenamestublen+2] = 'o';
   probdata->filenamestub[probdata->filenamestublen+3] = 'l';
   probdata->filenamestub[probdata->filenamestublen+4] = '\0';

   FILE* solfile = fopen(probdata->filenamestub, "w");
   if( solfile == NULL )
   {
      SCIPerrorMessage("could not open file <%s> for writing\n", probdata->filenamestub);
      probdata->filenamestub[probdata->filenamestublen] = '\0';

      return SCIP_WRITEERROR;
   }
   probdata->filenamestub[probdata->filenamestublen] = '\0';

   // see ampl/mp:sol.h:WriteSolFile() (seems buggy, https://github.com/ampl/mp/issues/135) and asl/writesol.c for solution file format
   SCIP_CALL( SCIPprintStatus(scip, solfile) );
   SCIPinfoMessage(scip, solfile, "\n\n");

   SCIPinfoMessage(scip, solfile, "Options\n%d\n", probdata->namplopts);
   for( int i = 0; i < probdata->namplopts; ++i )
      SCIPinfoMessage(scip, solfile, "%d\n", probdata->amplopts[i]);

   bool haveprimal = SCIPgetBestSol(scip) != NULL;
   bool havedual = probdata->islp && SCIPgetStage(scip) == SCIP_STAGE_SOLVED && !SCIPhasPerformedPresolve(scip);

   SCIPinfoMessage(scip, solfile, "%d\n%d\n", probdata->nconss, havedual ? probdata->nconss : 0);
   SCIPinfoMessage(scip, solfile, "%d\n%d\n", probdata->nvars, haveprimal ? probdata->nvars : 0);

   SCIPdebug( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, TRUE); )

   if( havedual )
      for( int c = 0; c < probdata->nconss; ++c )
      {
         SCIP_CONS* transcons;
         SCIP_Real dualval;

         /* dual solution is created by LP solver and therefore only available for linear constraints */
         SCIP_CALL( SCIPgetTransformedCons(scip, probdata->conss[c], &transcons) );
         assert(transcons == NULL || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(transcons)), "linear") == 0);

         if( transcons == NULL )
            dualval = 0.0;
         else if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
            dualval = SCIPgetDualsolLinear(scip, transcons);
         else
            dualval = -SCIPgetDualsolLinear(scip, transcons);
         assert(dualval != SCIP_INVALID);

         SCIPinfoMessage(scip, solfile, "%.17g\n", dualval);
      }

   if( haveprimal )
      for( int i = 0; i < probdata->nvars; ++i )
         SCIPinfoMessage(scip, solfile, "%.17g\n", SCIPgetSolVal(scip, SCIPgetBestSol(scip), probdata->vars[i]));

   /* AMPL solve status codes are at http://www.ampl.com/NEW/statuses.html
    *     number   string       interpretation
    *    0 -  99   solved       optimal solution found
    *  100 - 199   solved?      optimal solution indicated, but error likely
    *  200 - 299   infeasible   constraints cannot be satisfied
    *  300 - 399   unbounded    objective can be improved without limit
    *  400 - 499   limit        stopped by a limit that you set (such as on iterations)
    *  500 - 599   failure      stopped by an error condition in the solver routines
    */
   int solve_result_num;
   switch( SCIPgetStatus(scip) )
   {
      case SCIP_STATUS_UNKNOWN:
         solve_result_num = 500;
         break;
      case SCIP_STATUS_USERINTERRUPT:
         solve_result_num = 450;
         break;
      case SCIP_STATUS_NODELIMIT:
         solve_result_num = 400;
         break;
      case SCIP_STATUS_TOTALNODELIMIT:
         solve_result_num = 401;
         break;
      case SCIP_STATUS_STALLNODELIMIT:
         solve_result_num = 402;
         break;
      case SCIP_STATUS_TIMELIMIT:
         solve_result_num = 403;
         break;
      case SCIP_STATUS_MEMLIMIT:
         solve_result_num = 404;
         break;
      case SCIP_STATUS_GAPLIMIT:
         solve_result_num = 405;
         break;
      case SCIP_STATUS_SOLLIMIT:
         solve_result_num = 406;
         break;
      case SCIP_STATUS_BESTSOLLIMIT:
         solve_result_num = 407;
         break;
      case SCIP_STATUS_OPTIMAL:
         solve_result_num = 0;
         break;
      case SCIP_STATUS_INFEASIBLE:
         solve_result_num = 200;
         break;
      case SCIP_STATUS_UNBOUNDED:
         solve_result_num = 300;
         break;
      case SCIP_STATUS_INFORUNBD:
         solve_result_num = 299;
         break;
      default:
         /* solve_result_num = 500; */
         SCIPerrorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
         (void) fclose(solfile);
         return SCIP_INVALIDDATA;
   }
   SCIPinfoMessage(scip, solfile, "objno 0 %d\n", solve_result_num);

   if( fclose(solfile) != 0 )
   {
      SCIPerrorMessage("could not close solution file after writing\n");
      return SCIP_WRITEERROR;
   }

   return SCIP_OKAY;
}
