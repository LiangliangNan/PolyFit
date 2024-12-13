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

/**@file   reader_zpl.c
 * @ingroup DEFPLUGINS_READER
 * @brief  ZIMPL model file reader
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/reader_zpl.h"

#ifdef SCIP_WITH_ZIMPL

#include <unistd.h>
#include <stdbool.h>
#include <string.h>

#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/pub_misc.h"
#include "scip/pub_nlp.h"
#include "scip/pub_reader.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/cons_nonlinear.h"
#include "scip/struct_misc.h"
#include "scip/expr_pow.h"
#include "scip/expr_log.h"
#include "scip/expr_exp.h"
#include "scip/expr_abs.h"
#include "scip/expr_sum.h"
#include "scip/expr_trig.h"
#include "scip/expr_product.h"
#include "scip/pub_expr.h"
#include "scip/type_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

/* @Note: Due to dependencies we need the following order. */
/* include the ZIMPL headers necessary to define the LP and MINLP construction interface */
#include "zimpl/attribute.h"
#include "zimpl/ratlptypes.h"
#include "zimpl/lint.h"
#include "zimpl/mme.h"

#include "zimpl/numb.h"
#include "zimpl/bound.h"
#include "zimpl/mono.h"
#include "zimpl/term.h"

#include "zimpl/xlpglue.h"
#include "zimpl/zimpllib.h"

#ifdef __cplusplus
}
#endif

#define READER_NAME             "zplreader"
#define READER_DESC             "file reader for ZIMPL model files"
#define READER_EXTENSION        "zpl"

/*
 * LP construction interface of ZIMPL
 */

/* we only support ZIMPL with a version higher than 3.4.1 */
#if (ZIMPL_VERSION >= 341)

/* ZIMPL does not support user data in callbacks - we have to use static variables */
struct
SCIP_ReaderData
{
   SCIP*                 scip;               /**< scip data structure */
   SCIP_SOL*             sol;                /**< primal solution candidate */
   SCIP_Bool             valid;              /**< is the primal solution candidate valid */
   SCIP_Bool             branchpriowarning;  /**< store if the waring regarding fractional value for the branching
                                              *   priority was already posted */
   SCIP_Bool             initialconss;       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss;       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamiccols;        /**< should columns be added and removed dynamically to the LP? */
   SCIP_Bool             dynamicrows;        /**< should rows be added and removed dynamically to the LP? */
   SCIP_Bool             readerror;          /**< was a reading error be discovered */
   SCIP_RETCODE          retcode;            /**< store a none SCIP_OKAY return code if an error occurred */
};

/** create problem */
static
SCIP_RETCODE createProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   const char*           name                /**< name of the problem */
   )
{
   SCIP_Bool usestartsol;

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* check if are interested in the primal solution candidate */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/zplreader/usestartsol", &usestartsol) );

   if( usestartsol )
   {
      /* create primal solution */
      SCIP_CALL( SCIPcreateSol(scip, &readerdata->sol, NULL) );
      readerdata->valid = TRUE;
   }

   return SCIP_OKAY;
}

/** Allocate storage for the mathematical program instance generated by ZIMPL. xlp_alloc() is the first xlpglue routine
 *  that will be called by ZIMPL. The user_data pointer may hold an arbitray value.
 */
Lps* xlp_alloc(
   const char*           name,               /**< name of the problem */
   bool                  need_startval,      /**< does ZIMPL provides a primal solution candidate */
   void*                 user_data           /**< user data which was previously passed to ZIMPL */
   )
{  /*lint --e{715}*/
   SCIP* scip;
   SCIP_READERDATA* readerdata;

   readerdata = (SCIP_READERDATA*)user_data;
   assert(readerdata != NULL);
   assert(readerdata->retcode == SCIP_OKAY);
   assert(!readerdata->readerror);

   scip = readerdata->scip;
   assert(scip != NULL);

   readerdata->retcode = createProb(scip, readerdata, name);

   /* return the reader data pointer to receive it all other ZIMPL call backs */
   return (Lps*) readerdata;
}

/** free storage for mathematical program. xlp_free() is the last xlpglue routine that will be called by Zimpl */
void xlp_free(
   Lps*                  data                /**< pointer to reader data */
   )
{  /*lint --e{715}*/
   /* nothing to be done here */
}

/** does there already exists a constraint with the given name? */
bool xlp_conname_exists(
   const Lps*            data,               /**< pointer to reader data */
   const char*           name                /**< constraint name to check */
   )
{
   SCIP_READERDATA* readerdata;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   /* check if constraint with the given name already exists */
   return (SCIPfindCons(readerdata->scip, name) != NULL);
}

/** create a SCIP expression from a ZIMPL term
 *
 * Returns *expr == NULL if could not create expression due to unsupported ZIMPL functions.
 */
static
SCIP_RETCODE createExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   SCIP_EXPR**           expr,               /**< buffer to store expression */
   const Term*           term                /**< term to convert to expression */
)
{
   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(expr != NULL);
   assert(term != NULL);

   *expr = NULL;

   if( term_get_degree(term) == 2 )
   {
      int        nlinvars;
      int        nquadterms;
      SCIP_VAR** linvars;
      SCIP_VAR** quadvar1;
      SCIP_VAR** quadvar2;
      SCIP_Real* lincoefs;
      SCIP_Real* quadcoefs;
      Mono*      monom;
      int i;

      nlinvars   = 0;
      nquadterms = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &linvars,   term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvar1,  term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvar2,  term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs,  term_get_elements(term)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, term_get_elements(term)) );

      for( i = 0; i < term_get_elements(term); ++i )
      {
         monom = term_get_element(term, i);
         assert(!numb_equal(mono_get_coeff(monom), numb_zero()));
         assert(mono_get_degree(monom) <= 2);
         assert(mono_get_degree(monom) > 0);
         if (mono_get_degree(monom) == 1)
         {
            linvars [nlinvars] = (SCIP_VAR*)mono_get_var(monom, 0);
            lincoefs[nlinvars] = numb_todbl(mono_get_coeff(monom));
            ++nlinvars;
         }
         else
         {
            assert(mono_get_degree(monom) == 2);
            quadvar1 [nquadterms] = (SCIP_VAR*)mono_get_var(monom, 0);
            quadvar2 [nquadterms] = (SCIP_VAR*)mono_get_var(monom, 1);
            quadcoefs[nquadterms] = numb_todbl(mono_get_coeff(monom));
            ++nquadterms;
         }
      }

      SCIP_CALL( SCIPcreateExprQuadratic(scip, expr, nlinvars, linvars, lincoefs, nquadterms, quadvar1, quadvar2, quadcoefs, NULL, NULL) );

      SCIPfreeBufferArray(scip, &linvars);
      SCIPfreeBufferArray(scip, &quadvar1);
      SCIPfreeBufferArray(scip, &quadvar2);
      SCIPfreeBufferArray(scip, &lincoefs);
      SCIPfreeBufferArray(scip, &quadcoefs);
   }
   else
   {
      SCIP_VAR** polyvars;
      SCIP_Real* polyexps;
      SCIP_HASHMAP* varexpmap;
      SCIP_EXPR** monomials;
      int        nmonomials;
      int        monomialssize;
      SCIP_Real* coefs;
      Mono*      monomial;
      SCIP_EXPR* monomialexpr;
      SCIP_Bool created;
      int varpos;
      int i;
      int j;

      polyvars = NULL;
      polyexps = NULL;

      monomials = NULL;
      nmonomials = 0;
      monomialssize = 0;
      coefs = NULL;
      created = TRUE;

      SCIP_CALL( SCIPhashmapCreate(&varexpmap, SCIPblkmem(scip), SCIPcalcMemGrowSize(scip, 10)) );

      for( i = 0; i < term_get_elements(term); ++i )
      {
         monomial = term_get_element(term, i);
         assert(monomial != NULL);
         assert(!numb_equal(mono_get_coeff(monomial), numb_zero()));
         assert(mono_get_degree(monomial) > 0);

         /* allocate space in the monomials array */
         if( monomialssize == 0 )
         {
            monomialssize = SCIPcalcMemGrowSize(scip, 1);
            SCIP_CALL( SCIPallocBufferArray(scip, &monomials, monomialssize) );
            SCIP_CALL( SCIPallocBufferArray(scip, &coefs, monomialssize) );
         }
         else if( monomialssize < nmonomials + 1 )
         {
            monomialssize = SCIPcalcMemGrowSize(scip, nmonomials+1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &monomials, monomialssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, monomialssize) );
         }
         assert(monomials != NULL);
         assert(coefs != NULL);

         /* create SCIP monomial expression */
         for( j = 0; j < mono_get_degree(monomial); ++j )
         {
            SCIP_Real exponent;

            exponent = SCIPhashmapGetImageReal(varexpmap, (void*)mono_get_var(monomial, j));
            exponent = exponent == SCIP_INVALID ? 1.0 : exponent + 1.0;

            SCIP_CALL( SCIPhashmapSetImageReal(varexpmap, (void*)mono_get_var(monomial, j), exponent) );
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &polyvars, SCIPhashmapGetNElements(varexpmap)) );
         SCIP_CALL( SCIPallocBufferArray(scip, &polyexps, SCIPhashmapGetNElements(varexpmap)) );

         varpos = 0;

         for( j = 0; j < SCIPhashmapGetNEntries(varexpmap); ++j )
         {
            SCIP_HASHMAPENTRY* entry;

            entry = SCIPhashmapGetEntry(varexpmap, j);
            if( entry == NULL )
               continue;

            polyvars[varpos] = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
            polyexps[varpos] = SCIPhashmapEntryGetImageReal(entry);
            ++varpos;
         }
         assert(varpos == SCIPhashmapGetNElements(varexpmap));
         SCIPhashmapRemoveAll(varexpmap);

         SCIP_CALL( SCIPcreateExprMonomial(scip, &monomialexpr, varpos, polyvars, polyexps, NULL, NULL) );

         SCIPfreeBufferArrayNull(scip, &polyexps);
         SCIPfreeBufferArrayNull(scip, &polyvars);

         /* add monomial to array, possibly with an extra function around it */
         if( mono_get_function(monomial) == MFUN_NONE )
         {
            monomials[nmonomials] = monomialexpr;
            coefs[nmonomials] = numb_todbl(mono_get_coeff(monomial));
         }
         else
         {
            SCIP_EXPR* cosexpr;
            SCIP_EXPR* prodchildren[2];

            coefs[nmonomials] = 1.0;

            /* nonlinear monomial with an extra function around it */
            switch( mono_get_function(monomial) )
            {
            case MFUN_SQRT:
               SCIP_CALL( SCIPcreateExprPow(scip, &monomials[nmonomials], monomialexpr, 0.5, NULL, NULL) );
               break;
            case MFUN_LOG:
               /* log10(x) = ln(x) / ln(10.0) */
               coefs[nmonomials] = 1.0 / log(10.0);
               SCIP_CALL( SCIPcreateExprLog(scip, &monomials[nmonomials], monomialexpr, NULL, NULL) );
               break;
            case MFUN_EXP:
               SCIP_CALL( SCIPcreateExprExp(scip, &monomials[nmonomials], monomialexpr, NULL, NULL) );
               break;
            case MFUN_LN:
               SCIP_CALL( SCIPcreateExprLog(scip, &monomials[nmonomials], monomialexpr, NULL, NULL) );
               break;
            case MFUN_SIN:
               SCIP_CALL( SCIPcreateExprSin(scip, &monomials[nmonomials], monomialexpr, NULL, NULL) );
               break;
            case MFUN_COS:
               SCIP_CALL( SCIPcreateExprCos(scip, &monomials[nmonomials], monomialexpr, NULL, NULL) );
               break;
            case MFUN_TAN:
               SCIP_CALL( SCIPcreateExprSin(scip, &prodchildren[0], monomialexpr, NULL, NULL) );
               SCIP_CALL( SCIPcreateExprCos(scip, &cosexpr, monomialexpr, NULL, NULL) );
               SCIP_CALL( SCIPcreateExprPow(scip, &prodchildren[1], cosexpr, -1.0, NULL, NULL) );
               SCIP_CALL( SCIPcreateExprProduct(scip, &monomials[nmonomials], 2, prodchildren, 1.0, NULL, NULL) );

               SCIP_CALL( SCIPreleaseExpr(scip, &prodchildren[1]) );
               SCIP_CALL( SCIPreleaseExpr(scip, &cosexpr) );
               SCIP_CALL( SCIPreleaseExpr(scip, &prodchildren[0]) );

               break;
            case MFUN_ABS:
               SCIP_CALL( SCIPcreateExprAbs(scip, &monomials[nmonomials], monomialexpr, NULL, NULL) );
               break;
            case MFUN_POW:
               SCIP_CALL( SCIPcreateExprPow(scip, &monomials[nmonomials], monomialexpr,
                     numb_todbl(mono_get_coeff(monomial)), NULL, NULL) );
               break;
            case MFUN_SGNPOW:
               SCIP_CALL( SCIPcreateExprSignpower(scip, &monomials[nmonomials], monomialexpr,
                     numb_todbl(mono_get_coeff(monomial)), NULL, NULL) );
               break;
            case MFUN_NONE:
            case MFUN_TRUE:
            case MFUN_FALSE:
               SCIPerrorMessage("ZIMPL function %d invalid here.\n", mono_get_function(monomial));
               created = FALSE;
               break;
            default:
               SCIPerrorMessage("ZIMPL function %d not supported\n", mono_get_function(monomial));
               created = FALSE;
               break;
            }  /*lint !e788*/

            SCIP_CALL( SCIPreleaseExpr(scip, &monomialexpr) );
         }

         ++nmonomials;

         if( !created )
            break;
      }

      if( created )
      {
         SCIP_CALL( SCIPcreateExprSum(scip, expr, nmonomials, monomials, coefs, 0.0, NULL, NULL) );
      }

      /* free memory */
      for( j = nmonomials - 1; j >= 0; --j )
      {
         if( monomials[j] != NULL )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &monomials[j]) );
         }
      }

      SCIPfreeBufferArrayNull(scip, &coefs);
      SCIPfreeBufferArrayNull(scip, &monomials);
      SCIPhashmapFree(&varexpmap);
   }

   return SCIP_OKAY;
}

/** method creates a constraint and is called directly from ZIMPL
 *
 *  @note this method is used by ZIMPL beginning from version 3.00
 */
static
SCIP_RETCODE addConsTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   const char*           name,               /**< constraint name */
   ConType               type,               /**< constraint type (LHS, RHS, EQUAL, RANGE, etc) */
   const Numb*           lhs,                /**< left hand side */
   const Numb*           rhs,                /**< right hand side */
   unsigned int          flags,              /**< special constraint flags, see ratlptypes.h */
   const Term*           term,               /**< term to use */
   SCIP_Bool*            created             /**< pointer to store if a constraint was created */
   )
{
   SCIP_CONS* cons;
   SCIP_Real sciplhs;
   SCIP_Real sciprhs;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool usercut;
   SCIP_Bool lazycut;
   int i;

   switch( type )
   {
   case CON_FREE:
      sciplhs = -SCIPinfinity(scip);
      sciprhs = SCIPinfinity(scip);
      break;
   case CON_LHS:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = SCIPinfinity(scip);
      break;
   case CON_RHS:
      sciplhs = -SCIPinfinity(scip);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_RANGE:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      break;
   case CON_EQUAL:
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      assert(sciplhs == sciprhs);  /*lint !e777*/
      break;
   default:
      SCIPwarningMessage(scip, "invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
      sciplhs = (SCIP_Real)numb_todbl(lhs);
      sciprhs = (SCIP_Real)numb_todbl(rhs);
      readerdata->readerror = TRUE;
      break;
   }

   cons = NULL;

   /* default values */
   initial = readerdata->initialconss;
   separate = TRUE;
   propagate = TRUE;
   enforce = TRUE;
   check = TRUE;
   local = FALSE;
   modifiable = FALSE;

   usercut = (flags & LP_FLAG_CON_SEPAR) != 0;
   lazycut = (flags & LP_FLAG_CON_CHECK) != 0;

   /* evaluate constraint flags */
   if( usercut && lazycut )
   {
      initial = FALSE;
      separate = TRUE;
      check = TRUE;
   }
   else if( usercut )
   {
      initial = FALSE;
      separate = TRUE;
      check = FALSE;
   }
   else if( lazycut )
   {
      initial = FALSE;
      separate = FALSE;
      check = TRUE;
   }

   if( term_is_linear(term) )
   {
      /* if the constraint gives an indicator constraint */
      if ( flags & LP_FLAG_CON_INDIC )
      {
         bool lhsIndCons = FALSE;  /* generate lhs form for indicator constraints */
         bool rhsIndCons = FALSE;  /* generate rhs form for indicator constraints */

         /* currently indicator constraints can only handle "<=" constraints */
         switch( type )
         {
         case CON_LHS:
            lhsIndCons = TRUE;
            break;
         case CON_RHS:
            rhsIndCons = TRUE;
            break;
         case CON_RANGE:
         case CON_EQUAL:
            lhsIndCons = TRUE;
            rhsIndCons = TRUE;
            break;
         case CON_FREE:
            /*lint -fallthrough*/
         default:
            SCIPerrorMessage("invalid constraint type <%d> in ZIMPL callback xlp_addcon()\n", type);
            readerdata->readerror = TRUE;
            break;
         }

         /* insert lhs form of indicator */
         if ( lhsIndCons )
         {
            SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, NULL, 0, NULL, NULL, -sciplhs,
                  initial, separate, enforce, check, propagate, local, readerdata->dynamicconss, readerdata->dynamicrows, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );

            for( i = 0; i < term_get_elements(term); i++ )
            {
               SCIP_VAR* scipvar;
               SCIP_Real scipval;
               const Mono* mono = term_get_element(term, i);
               MFun mfun;

               scipvar = (SCIP_VAR*)mono_get_var(mono, 0);

               /* check whether variable is the binary variable */
               mfun = mono_get_function(mono);
               if (mfun == MFUN_TRUE || mfun == MFUN_FALSE)
               {
                  scipvar = (SCIP_VAR*)mono_get_var(mono, 0);
                  SCIP_CALL( SCIPsetBinaryVarIndicator(scip, cons, scipvar) );
               }
               else
               {
                  assert(!numb_equal(mono_get_coeff(mono), numb_zero()));
                  assert(mono_is_linear(mono));

                  scipval = -numb_todbl(mono_get_coeff(mono));
                  SCIP_CALL( SCIPaddVarIndicator(scip, cons, scipvar, scipval) );
               }
            }

            (*created) = TRUE;
         }

         /* insert rhs form of indicator */
         if ( rhsIndCons )
         {
            SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, NULL, 0, NULL, NULL, sciprhs,
                  initial, separate, enforce, check, propagate, local, readerdata->dynamicconss, readerdata->dynamicrows, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );

            for( i = 0; i < term_get_elements(term); i++ )
            {
               SCIP_VAR* scipvar;
               SCIP_Real scipval;
               const Mono* mono = term_get_element(term, i);
               MFun mfun;

               scipvar = (SCIP_VAR*)mono_get_var(mono, 0);

               /* check whether variable is the binary variable */
               mfun = mono_get_function(mono);
               if (mfun == MFUN_TRUE || mfun == MFUN_FALSE)
               {
                  scipvar = (SCIP_VAR*)mono_get_var(mono, 0);
                  SCIP_CALL( SCIPsetBinaryVarIndicator(scip, cons, scipvar) );
               }
               else
               {
                  assert(!numb_equal(mono_get_coeff(mono), numb_zero()));
                  assert(mono_is_linear(mono));

                  scipval = numb_todbl(mono_get_coeff(mono));
                  SCIP_CALL( SCIPaddVarIndicator(scip, cons, scipvar, scipval) );
               }
            }

            (*created) = TRUE;
         }
      }
      else
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 0, NULL, NULL, sciplhs, sciprhs,
               initial, separate, enforce, check, propagate, local, modifiable, readerdata->dynamicconss, readerdata->dynamicrows, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );

         for( i = 0; i < term_get_elements(term); i++ )
         {
            SCIP_VAR* scipvar;
            SCIP_Real scipval;

            assert(!numb_equal(mono_get_coeff(term_get_element(term, i)), numb_zero()));
            assert(mono_is_linear(term_get_element(term, i)));

            scipvar = (SCIP_VAR*)mono_get_var(term_get_element(term, i), 0);
            scipval = numb_todbl(mono_get_coeff(term_get_element(term, i)));

            SCIP_CALL( SCIPaddCoefLinear(scip, cons, scipvar, scipval) );
         }

         (*created) = TRUE;
      }
   }
   else
   {
      SCIP_EXPR* expr;

      /* convert term into expression */
      SCIP_CALL( createExpr(scip, readerdata, &expr, term) );

      if( expr == NULL )
      {
         /* ZIMPL term could not be represented as SCIP expression */
         (*created) = FALSE;
      }
      else
      {
         /* create constraint with expression */
         SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, name, expr, sciplhs, sciprhs,
            initial, separate, enforce, check, propagate, local, modifiable, readerdata->dynamicconss, readerdata->dynamicrows) );
         SCIP_CALL( SCIPaddCons(scip, cons) );

         SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

         (*created) = TRUE;
      }
   }

   if( cons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/** method adds objective term and is called directly from ZIMPL
 *
 *  @note this method is used by ZIMPL beginning from version 3.4.1
 */
static
SCIP_RETCODE addObjTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   const Term*           term                /**< term to use */
   )
{
   SCIP_Real objoffset;

   if( term_is_linear(term) )
   {
      int i;
      for( i = 0; i < term_get_elements(term); i++ )
      {
         SCIP_VAR* scipvar;
         SCIP_Real scipval;

         assert(!numb_equal(mono_get_coeff(term_get_element(term, i)), numb_zero()));
         assert(mono_is_linear(term_get_element(term, i)));

         scipvar = (SCIP_VAR*)mono_get_var(term_get_element(term, i), 0);
         scipval = numb_todbl(mono_get_coeff(term_get_element(term, i)));

         SCIP_CALL( SCIPaddVarObj(scip, scipvar, scipval) );
      }
   }
   else
   {
      /* create variable objvar, add 1*objvar to objective, and add constraint term - objvar = 0 */
      SCIP_EXPR* expr;
      SCIP_CONS* cons;
      SCIP_VAR* objvar;

      SCIP_CALL( createExpr(scip, readerdata, &expr, term) );

      if( expr == NULL )
      {
         SCIPerrorMessage("Could not convert ZIMPL objective term into SCIP expression due to unsupported ZIMPL function.\n");
         return SCIP_READERROR;
      }

      SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, "obj", expr,
         SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
         SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0,
         readerdata->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, readerdata->dynamicconss, FALSE) );

      SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, objvar, -1.0) );

      SCIP_CALL( SCIPaddVar(scip, objvar) );
      SCIP_CALL( SCIPaddCons(scip, cons) );

      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseVar(scip, &objvar) );
   }

   objoffset = numb_todbl(term_get_constant(term));
   SCIP_CALL( SCIPaddOrigObjoffset(scip, objoffset) );

   return SCIP_OKAY;
}

/** method creates a constraint and is called directly from ZIMPL
 *
 *  @note this method is used by ZIMPL beginning from version 3.00
 */
bool xlp_addcon_term(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name,               /**< constraint name */
   ConType               type,               /**< constraint type (LHS, RHS, EQUAL, RANGE, etc) */
   const Numb*           lhs,                /**< left hand side */
   const Numb*           rhs,                /**< right hand side */
   unsigned int          flags,              /**< special constraint flags, see ratlptypes.h */
   const Term*           term                /**< term to use */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_Bool created = FALSE;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   if( readerdata->retcode != SCIP_OKAY || readerdata->readerror )
      return TRUE;

   readerdata->retcode = addConsTerm(scip, readerdata, name, type, lhs, rhs, flags, term, &created);

   return !created;
}

/** adde variable */
static
SCIP_RETCODE addVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   const char*           name,               /**< variable name */
   VarClass              usevarclass,        /**< variable type */
   const Bound*          lower,              /**< lower bound */
   const Bound*          upper,              /**< upper bound */
   const Numb*           priority,           /**< branching priority */
   const Numb*           startval,           /**< start value for the variable within in the start solution */
   Var**                 zplvar              /**< pointer to store the created variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_VARTYPE vartype;
   SCIP_Bool initial;
   SCIP_Bool removable;
   int branchpriority;

   switch( bound_get_type(lower) )
   {
   case BOUND_VALUE:
      lb = (SCIP_Real)numb_todbl(bound_get_value(lower));
      break;
   case BOUND_INFTY:
      lb = SCIPinfinity(scip);
      break;
   case BOUND_MINUS_INFTY:
      lb = -SCIPinfinity(scip);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid lower bound type <%d> in ZIMPL reader\n", bound_get_type(lower));
      lb = 0.0;
      break;
   }

   switch( bound_get_type(upper) )
   {
   case BOUND_VALUE:
      ub = (SCIP_Real)numb_todbl(bound_get_value(upper));
      break;
   case BOUND_INFTY:
      ub = SCIPinfinity(scip);
      break;
   case BOUND_MINUS_INFTY:
      ub = -SCIPinfinity(scip);
      break;
   case BOUND_ERROR:
   default:
      SCIPerrorMessage("invalid upper bound type <%d> in ZIMPL reader\n", bound_get_type(upper));
      ub = 0.0;
      break;
   }

   switch( usevarclass )
   {
   case VAR_CON:
      vartype = SCIP_VARTYPE_CONTINUOUS;
      break;
   case VAR_INT:
      vartype = SCIP_VARTYPE_INTEGER;
      break;
   case VAR_IMP:
      vartype = SCIP_VARTYPE_IMPLINT;
      break;
   default:
      SCIPwarningMessage(scip, "invalid variable class <%d> in ZIMPL callback xlp_addvar()\n", usevarclass);
      vartype = SCIP_VARTYPE_CONTINUOUS;
      readerdata->readerror = TRUE;
      break;
   }
   initial = !(readerdata->dynamiccols);
   removable = readerdata->dynamiccols;

   /* create variable */
   SCIP_CALL( SCIPcreateVar(scip, &var, name, lb, ub, 0.0, vartype, initial, removable, NULL, NULL, NULL, NULL, NULL) );

   /* add variable to the problem; we are releasing the variable later */
   SCIP_CALL( SCIPaddVar(scip, var) );

   if( !numb_equal(priority, numb_unknown()) )
   {
      if( numb_is_int(priority) )
	 branchpriority = numb_toint(priority);
      else
      {
	 if( !readerdata->branchpriowarning )
	 {
	    SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,
	       "ZIMPL reader: fractional branching priorities in input - rounding down to integer values\n");
	    readerdata->branchpriowarning = TRUE;
	 }
	 branchpriority = (int)numb_todbl(priority);
      }

      /* change the branching priority of the variable */
      SCIP_CALL( SCIPchgVarBranchPriority(scip, var, branchpriority) );
   }

   /* check if we are willing to except a primal solution candidate */
   if( readerdata->valid )
   {
      /* if the number is unknown we have no valid primal solution candidate */
      if( numb_equal(startval, numb_unknown()) )
      {
         SCIPdebugMsg(scip, "primal solution candidate contains an unknown value for variable <%s>(%g)\n",
            SCIPvarGetName(var), (SCIP_Real)numb_todbl(startval));
         readerdata->valid = FALSE;
      }
      else
      {
         assert(readerdata->sol != NULL);
         SCIPdebugMsg(scip, "change solution solution <%p>: <%s> = <%g>\n",
            (void*)readerdata->sol, SCIPvarGetName(var), (SCIP_Real)numb_todbl(startval));

         /* set value within the primal solution candidate */
         SCIP_CALL( SCIPsetSolVal(scip, readerdata->sol, var, (SCIP_Real)numb_todbl(startval)) );
      }
   }

   /* copy the variable pointer before we release the variable */
   (*zplvar) = (Var*)var;

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** method adds a variable; is called directly by ZIMPL */
Var* xlp_addvar(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name,               /**< variable name */
   VarClass              usevarclass,        /**< variable type */
   const Bound*          lower,              /**< lower bound */
   const Bound*          upper,              /**< upper bound */
   const Numb*           priority,           /**< branching priority */
   const Numb*           startval            /**< start value for the variable within in the start solution */
   )
{  /*lint --e{715}*/
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   Var* zplvar;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   zplvar = NULL;

   if( readerdata->retcode != SCIP_OKAY || readerdata->readerror )
      return NULL;

   readerdata->retcode = addVar(scip, readerdata, name, usevarclass, lower, upper, priority, startval, &zplvar);

   return zplvar;
}

/** add a SOS constraint. Add a given a Zimpl term as an SOS constraint to the mathematical program */
static
SCIP_RETCODE addSOS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   const char*           name,               /**< constraint name */
   SosType               type,               /**< SOS type */
   const Term*           term                /**< terms indicating sos */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   int i;

   switch( type )
   {
   case SOS_TYPE1:
      separate = TRUE;
      enforce = TRUE;
      check = enforce;
      propagate = TRUE;
      local = FALSE;

      SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL,
            readerdata->initialconss, separate, enforce, check, propagate, local, readerdata->dynamicconss, readerdata->dynamicrows, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );

      for( i = 0; i < term_get_elements(term); i++ )
      {
         SCIP_VAR* var;
         SCIP_Real weight;

         assert( mono_is_linear(term_get_element(term, i)) );

         var = (SCIP_VAR*) mono_get_var(term_get_element(term, i), 0);
         weight = numb_todbl(mono_get_coeff(term_get_element(term, i)));

         SCIP_CALL( SCIPaddVarSOS1(scip, cons, var, weight) );
      }
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      break;
   case SOS_TYPE2:
      separate = TRUE;
      enforce = TRUE;
      check = enforce;
      propagate = TRUE;
      local = FALSE;

      SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL,
            readerdata->initialconss, separate, enforce, check, propagate, local, readerdata->dynamicconss, readerdata->dynamicrows, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      for( i = 0; i < term_get_elements(term); i++ )
      {
         SCIP_VAR* var;
         SCIP_Real weight;

         assert( mono_is_linear(term_get_element(term, i)) );

         var = (SCIP_VAR*) mono_get_var(term_get_element(term, i), 0);
         weight = numb_todbl(mono_get_coeff(term_get_element(term, i)));

         SCIP_CALL( SCIPaddVarSOS2(scip, cons, var, weight) );
      }
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      break;
   case SOS_ERR:
      /*lint -fallthrough*/
   default:
      SCIPerrorMessage("invalid SOS type <%d> in ZIMPL callback xlp_addsos_term()\n", type);
      readerdata->readerror = TRUE;
      break;
   }

   return SCIP_OKAY;
}

/** add a SOS constraint. Add a given a Zimpl term as an SOS constraint to the mathematical program */
int xlp_addsos_term(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name,               /**< constraint name */
   SosType               type,               /**< SOS type */
   const Numb*           priority,           /**< priority */
   const Term*           term                /**< terms indicating sos */
   )
{
   /*lint --e{715}*/
   SCIP* scip;
   SCIP_READERDATA* readerdata;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   if( readerdata->retcode != SCIP_OKAY || readerdata->readerror )
      return TRUE;

   readerdata->retcode = addSOS(scip, readerdata, name, type, term);

   return 0;
}

/** returns the variable name */
const char* xlp_getvarname(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
#ifndef NDEBUG
   SCIP* scip;
   SCIP_READERDATA* readerdata;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);
#endif

   return SCIPvarGetName((SCIP_VAR*)var);
}

/** return variable type */
VarClass xlp_getclass(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_VAR* scipvar;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scipvar = (SCIP_VAR*)var;
   switch( SCIPvarGetType(scipvar) )
   {
   case SCIP_VARTYPE_BINARY:
   case SCIP_VARTYPE_INTEGER:
      return VAR_INT;
   case SCIP_VARTYPE_IMPLINT:
      return VAR_IMP;
   case SCIP_VARTYPE_CONTINUOUS:
      return VAR_CON;
   default:
      SCIPerrorMessage("invalid SCIP variable type <%d> in ZIMPL callback xlp_getclass()\n", SCIPvarGetType(scipvar));
      readerdata->readerror = TRUE;
      break;
   }

   return VAR_CON;
}

/** returns lower bound */
Bound* xlp_getlower(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_VAR* scipvar;
   SCIP_Real lb;
   char s[SCIP_MAXSTRLEN];
   BoundType boundtype;
   Numb* numb;
   Bound* bound;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   scipvar = (SCIP_VAR*)var;
   assert(scipvar != NULL);

   /* collect lower bound */
   lb = SCIPvarGetLbGlobal(scipvar);
   numb = NULL;

   /* check if lower bound is infinity */
   if( SCIPisInfinity(scip, -lb) )
      boundtype = BOUND_MINUS_INFTY;
   else if( SCIPisInfinity(scip, lb) )
      boundtype = BOUND_INFTY;
   else
   {
      boundtype = BOUND_VALUE;

      /* create double form string */
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", lb);
      numb = numb_new_ascii(s);
   }

   /* create bound */
   bound = bound_new(boundtype, numb);

   if( numb != NULL )
      numb_free(numb);

   return bound;
}

/** returns upper bound */
Bound* xlp_getupper(
   const Lps*            data,               /**< pointer to reader data */
   const Var*            var                 /**< variable */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_VAR* scipvar;
   SCIP_Real ub;
   char s[SCIP_MAXSTRLEN];
   BoundType boundtype;
   Numb* numb;
   Bound* bound;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   scipvar = (SCIP_VAR*)var;
   assert(scipvar != NULL);

   /* collect upper bound */
   ub = SCIPvarGetUbGlobal(scipvar);
   numb = NULL;

   /* check if upper bound is infinity */
   if( SCIPisInfinity(scip, -ub) )
      boundtype = BOUND_MINUS_INFTY;
   else if( SCIPisInfinity(scip, ub) )
      boundtype = BOUND_INFTY;
   else
   {
      boundtype = BOUND_VALUE;
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%.20f", ub);
      numb = numb_new_ascii(s);
   }

   /* create ZIMPL bound */
   bound = bound_new(boundtype, numb);

   if (numb != NULL)
      numb_free(numb);

   return bound;
}

/** Set the name and direction of the objective function, i.e. minimization or maximization
 *  Coefficents of the objective function will be set to all zero.
 */
bool xlp_setobj(
   Lps*                  data,               /**< pointer to reader data */
   const char*           name,               /**< name of the objective function */
   bool                  minimize            /**< True if the problem should be minimized, False if it should be maximized  */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;
   SCIP_OBJSENSE objsense;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   if( readerdata->retcode != SCIP_OKAY || readerdata->readerror )
      return FALSE;

   objsense = (minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE);
   readerdata->retcode = SCIPsetObjsense(scip, objsense);

   return FALSE;
}

/** adds objective function */
void xlp_addtoobj(
   Lps*                  data,               /**< pointer to reader data */
   const Term*           term                /**< objective term */
   )
{
   SCIP* scip;
   SCIP_READERDATA* readerdata;

   readerdata = (SCIP_READERDATA*)data;
   assert(readerdata != NULL);

   scip = readerdata->scip;
   assert(scip != NULL);

   if( readerdata->retcode != SCIP_OKAY || readerdata->readerror )
      return;

   readerdata->retcode = addObjTerm(scip, readerdata, term);
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyZpl)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderZpl(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadZpl)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;
   SCIP_RETCODE retcode;
   char oldpath[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   char compextension[SCIP_MAXSTRLEN];
   char namewithoutpath[SCIP_MAXSTRLEN];
   char* path;
   char* name;
   char* extension;
   char* compression;
   char* paramstr;

   SCIP_Bool changedir;
   int i;

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/zplreader/changedir", &changedir) );

   path = NULL;
   oldpath[0] = '\0';
   if( changedir )
   {
      /* change to the directory of the ZIMPL file, s.t. paths of data files read by the ZIMPL model are relative to
       * the location of the ZIMPL file
       */
      (void)SCIPstrncpy(buffer, filename, SCIP_MAXSTRLEN);
      SCIPsplitFilename(buffer, &path, &name, &extension, &compression);
      if( compression != NULL )
         (void) SCIPsnprintf(compextension, SCIP_MAXSTRLEN, ".%s", compression);
      else
         *compextension = '\0';
      (void) SCIPsnprintf(namewithoutpath, SCIP_MAXSTRLEN, "%s.%s%s", name, extension, compextension);
      if( (char*)getcwd(oldpath, SCIP_MAXSTRLEN) == NULL )
      {
         SCIPerrorMessage("error getting the current path\n");
         return SCIP_READERROR;
      }
      if( path != NULL )
      {
         if( chdir(path) != 0 )
         {
            SCIPerrorMessage("error changing to directory <%s>\n", path);
            return SCIP_NOFILE;
         }
      }
      filename = namewithoutpath;
   }

   /* get current path for output */
   if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL )
   {
      char currentpath[SCIP_MAXSTRLEN];
      if( (char*)getcwd(currentpath, SCIP_MAXSTRLEN) == NULL )
      {
         SCIPerrorMessage("error getting the current path\n");
         return SCIP_READERROR;
      }
      /* an extra blank line should be printed separately since the buffer message handler only handle up to one line
       *  correctly */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "base directory for ZIMPL parsing: <%s>\n", currentpath);
      /* an extra blank line should be printed separately since the buffer message handler only handle up to one line
       *  correctly */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
   }

   /* allocate storage */
   SCIP_CALL( SCIPallocBuffer(scip, &readerdata) );

   readerdata->scip = scip;
   readerdata->sol = NULL;
   readerdata->valid = FALSE;
   readerdata->branchpriowarning = FALSE;
   readerdata->readerror = FALSE;
   readerdata->retcode = SCIP_OKAY;
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &(readerdata->initialconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicconss", &(readerdata->dynamicconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &(readerdata->dynamiccols)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &(readerdata->dynamicrows)) );

   /* get the parameter string */
   SCIP_CALL( SCIPgetStringParam(scip, "reading/zplreader/parameters", &paramstr) );
   if( strcmp(paramstr, "-") == 0 )
   {
      /* call ZIMPL parser without arguments */
      if( !zpl_read(filename, TRUE, (void*)readerdata) )
         readerdata->readerror = TRUE;
      else
      {
         /* evaluate retcode */
         if ( readerdata->retcode != SCIP_OKAY )
         {
            SCIPfreeBuffer(scip, &readerdata);
            return readerdata->retcode;
         }
      }
   }
   else
   {
      char dummy[2] = "x";
      char** argv;
      int argc;
      int p;
      int len;

      len = (int) strlen(paramstr);
      SCIP_CALL( SCIPallocBufferArray(scip, &argv, len+1) );
      argv[0] = dummy; /* argument 0 is irrelevant */
      argc = 1;
      p = 0;
      while( p < len )
      {
         int arglen;

         /* process next argument */
         SCIP_CALL( SCIPallocBufferArray(scip, &argv[argc], len+1) );  /*lint !e866*/
         arglen = 0;

         /* skip spaces */
         while( p < len && paramstr[p] == ' ' )
            p++;

         /* process characters */
         while( p < len && paramstr[p] != ' ' )
         {
            switch( paramstr[p] )
            {
            case '"':
               p++;
               /* read characters as they are until the next " */
               while( p < len && paramstr[p] != '"' )
               {
                  argv[argc][arglen] = paramstr[p];
                  arglen++;
                  p++;
               }
               p++; /* skip final " */
               break;
            case '\\':
               /* read next character as it is */
               p++;
               argv[argc][arglen] = paramstr[p];
               arglen++;
               p++;
               break;
            default:
               argv[argc][arglen] = paramstr[p];
               arglen++;
               p++;
               break;
            }
         }
         argv[argc][arglen] = '\0';

         /* check for empty argument */
         if( arglen == 0 )
         {
            SCIPfreeBufferArray(scip, &argv[argc]);
         }
         else
            argc++;
      }

      /* append file name as last argument */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &argv[argc], filename, (int) strlen(filename)+1) );  /*lint !e866*/
      argc++;

      /* display parsed arguments */
      if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "ZIMPL arguments:\n");
         for( i = 1; i < argc; ++i )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%d: <%s>\n", i, argv[i]);
         }
      }

      /* call ZIMPL parser with arguments */
      if( !zpl_read_with_args(argv, argc, TRUE, (void*)readerdata) )
         readerdata->readerror = TRUE;

      /* free argument memory */
      for( i = argc - 1; i >= 1; --i )
      {
         SCIPfreeBufferArray(scip, &argv[i]);
      }
      SCIPfreeBufferArray(scip, &argv);

      if ( readerdata->retcode != SCIP_OKAY )
      {
         SCIPfreeBuffer(scip, &readerdata);
         return readerdata->retcode;
      }
   }

   if( changedir )
   {
      /* change directory back to old path */
      if( path != NULL )
      {
         if( chdir(oldpath) != 0 )
         {
            SCIPwarningMessage(scip, "error changing back to directory <%s>\n", oldpath);
         }
      }
   }

   if( readerdata->valid )
   {
      SCIP_Bool stored;

      assert(readerdata->sol != NULL);

      stored = FALSE;

      /* add primal solution to solution candidate storage, frees the solution afterwards */
      SCIP_CALL( SCIPaddSolFree(scip, &readerdata->sol, &stored) );

      if( stored )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "ZIMPL starting solution candidate accepted\n");
      }
   }

   *result = SCIP_SUCCESS;

   /* evaluate if a reading error occurred */
   if( readerdata->readerror )
      retcode = SCIP_READERROR;
   else
      retcode = SCIP_OKAY;

   /* free primal solution candidate */
   if( readerdata->sol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &readerdata->sol) );
   }

   /* free reader data */
   SCIPfreeBuffer(scip, &readerdata);

   return retcode;
}


#endif
#endif


/*
 * reader specific interface methods
 */

/** includes the zpl file reader in SCIP */ /*lint --e{715}*/
SCIP_RETCODE SCIPincludeReaderZpl(
   SCIP*                 scip                /**< SCIP data structure */
   )
{ /*lint --e{715}*/
#ifdef SCIP_WITH_ZIMPL
#if (ZIMPL_VERSION >= 341)
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;
   char extcodename[SCIP_MAXSTRLEN];

   assert(scip != NULL);

   /* create zpl reader data */
   readerdata = NULL;
   reader = NULL;
   /* include zpl reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyZpl) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadZpl) );

   /* add zpl reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/zplreader/changedir", "should the current directory be changed to that of the ZIMPL file before parsing?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/zplreader/usestartsol", "should ZIMPL starting solutions be forwarded to SCIP?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddStringParam(scip,
         "reading/zplreader/parameters", "additional parameter string passed to the ZIMPL parser (or - for no additional parameters)",
         NULL, FALSE, "-", NULL, NULL) );

   (void) SCIPsnprintf(extcodename, SCIP_MAXSTRLEN, "ZIMPL %d.%d.%d", ZIMPL_VERSION/100, (ZIMPL_VERSION%100)/10, ZIMPL_VERSION%10); /*lint !e778*/
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, extcodename, "Zuse Institute Mathematical Programming Language developed by T. Koch (zimpl.zib.de)"));
#else
   assert(scip != NULL);

   SCIPwarningMessage(scip, "SCIP does only support ZIMPL 3.4.1 and higher. Please update your ZIMPL version %d.%d.%d\n",
      ZIMPL_VERSION/100, (ZIMPL_VERSION%100)/10, ZIMPL_VERSION%10);
#endif
#endif

   return SCIP_OKAY;
}
