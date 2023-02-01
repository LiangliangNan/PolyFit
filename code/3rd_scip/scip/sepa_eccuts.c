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

/**@file   sepa_eccuts.c
 * @brief  edge concave cut separator
 * @author Benjamin Müller
 */

/**@todo only count number of fixed variables in the edge concave terms */
/**@todo only add nonlinear row aggregations where at least ...% of the variables (bilinear terms?) are in edge concave
 * terms */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"
#include "scip/sepa_eccuts.h"
#include "scip/cons_xor.h"
#include "scip/nlp.h"
#include "tclique/tclique.h"

#define SEPA_NAME                            "eccuts"
#define SEPA_DESC                            "separator for edge-concave functions"
#define SEPA_PRIORITY            -13000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CLIQUE_MAXFIRSTNODEWEIGHT  1000 /**< maximum weight of branching nodes in level 0; 0 if not used for cliques
                                         *   with at least one fractional node) */
#define CLIQUE_MINWEIGHT              0 /**< lower bound for weight of generated cliques */
#define CLIQUE_MAXNTREENODES      10000 /**< maximal number of nodes of b&b tree */
#define CLIQUE_BACKTRACKFREQ      10000 /**< frequency to backtrack to first level of tree (0: no premature backtracking) */

#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MAXROUNDS            10 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT       250 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXDEPTH             -1 /**< maximal depth at which the separator is applied */
#define DEFAULT_MAXSEPACUTS          10 /**< maximal number of e.c. cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT      50 /**< maximal number of e.c. cuts separated per separation round in root node */
#define DEFAULT_CUTMAXRANGE        1e+7 /**< maximal coefficient range of a cut (maximal coefficient divided by minimal
                                         *   coefficient) in order to be added to LP relaxation */
#define DEFAULT_MINVIOLATION        0.3 /**< minimal violation of an e.c. cut to be separated */
#define DEFAULT_MINAGGRSIZE           3 /**< search for e.c. aggregation of at least this size (has to be >= 3) */
#define DEFAULT_MAXAGGRSIZE           4 /**< search for e.c. aggregation of at most this size (has to be >= minaggrsize) */
#define DEFAULT_MAXBILINTERMS       500 /**< maximum number of bilinear terms allowed to be in a quadratic constraint */
#define DEFAULT_MAXSTALLROUNDS        5 /**< maximum number of unsuccessful rounds in the e.c. aggregation search */

#define SUBSCIP_NODELIMIT         100LL /**< node limit to solve the sub-SCIP */

#define ADJUSTFACETTOL             1e-6 /**< adjust resulting facets in checkRikun() up to a violation of this value */
#define USEDUALSIMPLEX             TRUE /**< use dual or primal simplex algorithm? */

/** first values for 2^n */
static const int poweroftwo[] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192 };

/*
 * Data structures
 */

/** data to store a single edge-concave aggregations; an edge-concave aggregation of a quadratic constraint is a subset
 *  of nonconvex bilinear terms
 */
struct EcAggr
{
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   int                   varsize;            /**< size of vars array */

   SCIP_Real*            termcoefs;          /**< coefficients of bilinear terms */
   int*                  termvars1;          /**< index of the first variable of each bilinear term */
   int*                  termvars2;          /**< index of the second variable of each bilinear term*/
   int                   nterms;             /**< number of bilinear terms in the aggregation */
   int                   termsize;           /**< size of term{coefs,vars1,vars2} arrays */
};
typedef struct EcAggr SCIP_ECAGGR;

/** data to store all edge-concave aggregations and the remaining part of a nonlinear row of the form g(x) <= rhs */
struct NlrowAggr
{
   SCIP_NLROW*           nlrow;              /**< nonlinear row aggregation */
   SCIP_Bool             rhsaggr;            /**< consider nonlinear row aggregation for g(x) <= rhs (TRUE) or
                                              *   g(x) >= lhs (FALSE) */

   SCIP_ECAGGR**         ecaggr;             /**< array with all edge-concave aggregations */
   int                   necaggr;            /**< number of edge-concave aggregation */

   SCIP_VAR**            linvars;            /**< linear variables */
   SCIP_Real*            lincoefs;           /**< linear coefficients */
   int                   nlinvars;           /**< number of linear variables */

   SCIP_VAR**            quadvars;           /**< quadratic variables */
   int*                  quadvar2aggr;       /**< stores in which edge-concave aggregation the i-th quadratic variable
                                              *   is contained (< 0: in no edge-concave aggregation) */
   int                   nquadvars;          /**< number of quadratic variables */

   SCIP_VAR**            remtermvars1;       /**< first quadratic variable of remaining bilinear terms */
   SCIP_VAR**            remtermvars2;       /**< second quadratic variable of remaining bilinear terms */
   SCIP_Real*            remtermcoefs;       /**< coefficients for each remaining bilinear term */
   int                   nremterms;          /**< number of remaining bilinear terms */
   int                   remtermsize;        /**< size of remterm* arrays */

   SCIP_Real             rhs;                /**< rhs of the nonlinear row */
   SCIP_Real             constant;           /**< constant part of the nonlinear row */
};
typedef struct NlrowAggr SCIP_NLROWAGGR;

/** separator data */
struct SCIP_SepaData
{
   SCIP_NLROWAGGR**      nlrowaggrs;         /**< array containing all nonlinear row aggregations */
   int                   nnlrowaggrs;        /**< number of nonlinear row aggregations */
   int                   nlrowaggrssize;     /**< size of nlrowaggrs array */
   SCIP_Bool             searchedforaggr;    /**< flag if we already searched for nlrow aggregation candidates */
   int                   minaggrsize;        /**< only search for e.c. aggregations of at least this size (has to be >= 3) */
   int                   maxaggrsize;        /**< only search for e.c. aggregations of at most this size (has to be >= minaggrsize) */
   int                   maxecsize;          /**< largest edge concave aggregation size */
   int                   maxbilinterms;      /**< maximum number of bilinear terms allowed to be in a quadratic constraint */
   int                   maxstallrounds;     /**< maximum number of unsuccessful rounds in the e.c. aggregation search */

   SCIP_LPI*             lpi;                /**< LP interface to solve the LPs to compute the facets of the convex envelopes */
   int                   lpisize;            /**< maximum size of e.c. aggregations which can be handled by the LP interface */

   SCIP_Real             cutmaxrange;        /**< maximal coef range of a cut (maximal coefficient divided by minimal
                                              *   coefficient) in order to be added to LP relaxation */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Real             minviolation;       /**< minimal violation of an e.c. cut to be separated */

   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxdepth;           /**< maximal depth at which the separator is applied */
   int                   maxsepacuts;        /**< maximal number of e.c. cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of e.c. cuts separated per separation round in root node */

#ifdef SCIP_STATISTIC
   SCIP_Real             aggrsearchtime;     /**< total time spent for searching edge concave aggregations */
   int                   nlhsnlrowaggrs;     /**< number of found nonlinear row aggregations for SCIP_NLROWs of the form g(x) <= rhs */
   int                   nrhsnlrowaggrs;     /**< number of found nonlinear row aggregations for SCIP_NLROWs of the form g(x) >= lhs */
#endif
};


/*
 * Local methods
 */

/** creates and empty edge-concave aggregation (without bilinear terms) */
static
SCIP_RETCODE ecaggrCreateEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ECAGGR**         ecaggr,             /**< pointer to store the edge-concave aggregation */
   int                   nquadvars,          /**< number of quadratic variables */
   int                   nquadterms          /**< number of bilinear terms */
   )
{
   assert(scip != NULL);
   assert(ecaggr != NULL);
   assert(nquadvars > 0);
   assert(nquadterms >= nquadvars);

   SCIP_CALL( SCIPallocBlockMemory(scip, ecaggr) );

   (*ecaggr)->nvars = 0;
   (*ecaggr)->nterms = 0;
   (*ecaggr)->varsize = nquadvars;
   (*ecaggr)->termsize = nquadterms;

   /* allocate enough memory for the quadratic variables and bilinear terms */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*ecaggr)->vars, nquadvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*ecaggr)->termcoefs, nquadterms) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*ecaggr)->termvars1, nquadterms) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*ecaggr)->termvars2, nquadterms) );

   return SCIP_OKAY;
}

/** frees and edge-concave aggregation */
static
SCIP_RETCODE ecaggrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ECAGGR**         ecaggr              /**< pointer to store the edge-concave aggregation */
   )
{
   assert(scip != NULL);
   assert(ecaggr != NULL);

   SCIPfreeBlockMemoryArray(scip, &((*ecaggr)->termcoefs), (*ecaggr)->termsize);
   SCIPfreeBlockMemoryArray(scip, &((*ecaggr)->termvars1), (*ecaggr)->termsize);
   SCIPfreeBlockMemoryArray(scip, &((*ecaggr)->termvars2), (*ecaggr)->termsize);
   SCIPfreeBlockMemoryArray(scip, &((*ecaggr)->vars), (*ecaggr)->varsize);

   SCIPfreeBlockMemory(scip, ecaggr);
   *ecaggr = NULL;

   return SCIP_OKAY;
}

/** adds a quadratic variable to an edge-concave aggregation */
static
SCIP_RETCODE ecaggrAddQuadvar(
   SCIP_ECAGGR*          ecaggr,             /**< pointer to store the edge-concave aggregation */
   SCIP_VAR*             x                   /**< first variable */
   )
{
   ecaggr->vars[ ecaggr->nvars++ ] = x;
   return SCIP_OKAY;
}

/** adds a bilinear term to an edge-concave aggregation */
static
SCIP_RETCODE ecaggrAddBilinTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ECAGGR*          ecaggr,             /**< pointer to store the edge-concave aggregation */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             coef                /**< bilinear coefficient */
   )
{
   int idx1;
   int idx2;
   int i;

   assert(x != NULL);
   assert(y != NULL);
   assert(ecaggr->nterms + 1 <= ((ecaggr->nvars + 1) * ecaggr->nvars) / 2);
   assert(!SCIPisZero(scip, coef));

   idx1 = -1;
   idx2 = -1;

   /* search for the quadratic variables in the e.c. aggregation */
   for( i = 0; i < ecaggr->nvars && (idx1 == -1 || idx2 == -1); ++i )
   {
      if( ecaggr->vars[i] == x )
         idx1 = i;
      if( ecaggr->vars[i] == y )
         idx2 = i;
   }

   assert(idx1 != -1 && idx2 != -1);

   ecaggr->termcoefs[ ecaggr->nterms ] = coef;
   ecaggr->termvars1[ ecaggr->nterms ] = idx1;
   ecaggr->termvars2[ ecaggr->nterms ] = idx2;
   ++(ecaggr->nterms);

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** prints an edge-concave aggregation */
static
void ecaggrPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ECAGGR*          ecaggr              /**< pointer to store the edge-concave aggregation */
   )
{
   int i;

   assert(scip != NULL);
   assert(ecaggr != NULL);

   SCIPdebugMsg(scip, " nvars = %d nterms = %d\n", ecaggr->nvars, ecaggr->nterms);
   SCIPdebugMsg(scip, " vars: ");
   for( i = 0; i < ecaggr->nvars; ++i )
      SCIPdebugMsgPrint(scip, "%s ", SCIPvarGetName(ecaggr->vars[i]));
   SCIPdebugMsgPrint(scip, "\n");

   SCIPdebugMsg(scip, " terms: ");
   for( i = 0; i < ecaggr->nterms; ++i )
   {
      SCIP_VAR* x;
      SCIP_VAR* y;

      x = ecaggr->vars[ ecaggr->termvars1[i] ];
      y = ecaggr->vars[ ecaggr->termvars2[i] ];
      SCIPdebugMsgPrint(scip, "%e %s * %s  ", ecaggr->termcoefs[i], SCIPvarGetName(x), SCIPvarGetName(y) );
   }
   SCIPdebugMsgPrint(scip, "\n");
}
#endif

/** stores linear terms in a given nonlinear row aggregation */
static
SCIP_RETCODE nlrowaggrStoreLinearTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROWAGGR*       nlrowaggr,          /**< nonlinear row aggregation */
   SCIP_VAR**            linvars,            /**< linear variables */
   SCIP_Real*            lincoefs,           /**< linear coefficients */
   int                   nlinvars            /**< number of linear variables */
   )
{
   assert(scip != NULL);
   assert(nlrowaggr != NULL);
   assert(linvars != NULL || nlinvars == 0);
   assert(lincoefs != NULL || nlinvars == 0);
   assert(nlinvars >= 0);

   nlrowaggr->nlinvars = nlinvars;
   nlrowaggr->linvars = NULL;
   nlrowaggr->lincoefs = NULL;

   if( nlinvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlrowaggr->linvars, nlinvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlrowaggr->lincoefs, nlinvars) );
      BMScopyMemoryArray(nlrowaggr->linvars, linvars, nlinvars);
      BMScopyMemoryArray(nlrowaggr->lincoefs, lincoefs, nlinvars);
   }

   /* if we have a nlrow of the form g(x) >= lhs, multiply every coefficient by -1 */
   if( !nlrowaggr->rhsaggr )
   {
      int i;

      for( i = 0; i < nlrowaggr->nlinvars; ++i )
         nlrowaggr->lincoefs[i] *= -1.0;
   }

   return SCIP_OKAY;
}

/** stores quadratic variables in a given nonlinear row aggregation */
static
SCIP_RETCODE nlrowaggrStoreQuadraticVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROWAGGR*       nlrowaggr,          /**< nonlinear row aggregation */
   SCIP_VAR**            quadvars,           /**< quadratic variables */
   int                   nquadvars           /**< number of quadratic variables */
   )
{
   assert(scip != NULL);
   assert(nlrowaggr != NULL);
   assert(quadvars != NULL);
   assert(nquadvars > 0);

   nlrowaggr->nquadvars = nquadvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlrowaggr->quadvars, nquadvars) );
   BMScopyMemoryArray(nlrowaggr->quadvars, quadvars, nquadvars);

   return SCIP_OKAY;
}

/** adds a remaining bilinear term to a given nonlinear row aggregation */
static
SCIP_RETCODE nlrowaggrAddRemBilinTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROWAGGR*       nlrowaggr,          /**< nonlinear row aggregation */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             coef                /**< bilinear coefficient */
   )
{
   assert(nlrowaggr != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(coef != 0.0);

   nlrowaggr->remtermcoefs[ nlrowaggr->nremterms ] = coef;
   nlrowaggr->remtermvars1[ nlrowaggr->nremterms ] = x;
   nlrowaggr->remtermvars2[ nlrowaggr->nremterms ] = y;
   ++(nlrowaggr->nremterms);

   return SCIP_OKAY;
}

/** creates a nonlinear row aggregation */
static
SCIP_RETCODE nlrowaggrCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_NLROWAGGR**      nlrowaggr,          /**< pointer to store the nonlinear row aggregation */
   int*                  quadvar2aggr,       /**< mapping between quadratic variables and edge-concave aggregation
                                              *   stores a negative value if the quadratic variables does not belong
                                              *   to any aggregation */
   int                   nfound,             /**< number of edge-concave aggregations */
   SCIP_Bool             rhsaggr             /**< consider nonlinear row aggregation for g(x) <= rhs (TRUE) or
                                              *   lhs <= g(x) (FALSE) */
   )
{
   int* aggrnvars;  /* count the number of variables in each e.c. aggregations */
   int* aggrnterms; /* count the number of bilinear terms in each e.c. aggregations */
   int nquadelems;
   int nquadvars;
   int nremterms;
   int i;

   assert(scip != NULL);
   assert(nlrow != NULL);
   assert(nlrowaggr != NULL);
   assert(quadvar2aggr != NULL);
   assert(nfound > 0);

   nquadelems = SCIPnlrowGetNQuadElems(nlrow);
   nquadvars = SCIPnlrowGetNQuadVars(nlrow);
   nremterms = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &aggrnvars, nfound) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrnterms, nfound) );

   /* create an empty nonlinear row aggregation */
   SCIP_CALL( SCIPallocBlockMemory(scip, nlrowaggr) );
   (*nlrowaggr)->nlrow = nlrow;
   (*nlrowaggr)->rhsaggr = rhsaggr;
   (*nlrowaggr)->nquadvars = nquadvars;
   (*nlrowaggr)->rhs = rhsaggr ? SCIPnlrowGetRhs(nlrow) : -SCIPnlrowGetLhs(nlrow);
   (*nlrowaggr)->constant = rhsaggr ? SCIPnlrowGetConstant(nlrow) : -SCIPnlrowGetConstant(nlrow);

   (*nlrowaggr)->quadvars = NULL;
   (*nlrowaggr)->quadvar2aggr = NULL;
   (*nlrowaggr)->remtermcoefs = NULL;
   (*nlrowaggr)->remtermvars1 = NULL;
   (*nlrowaggr)->remtermvars2 = NULL;
   (*nlrowaggr)->nquadvars = 0;
   (*nlrowaggr)->nremterms = 0;

   /* copy quadvar2aggr array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlrowaggr)->quadvar2aggr, nquadvars) );
   BMScopyMemoryArray((*nlrowaggr)->quadvar2aggr, quadvar2aggr, nquadvars);

   /* store all linear terms */
   SCIP_CALL( nlrowaggrStoreLinearTerms(scip, *nlrowaggr, SCIPnlrowGetLinearVars(nlrow), SCIPnlrowGetLinearCoefs(nlrow),
         SCIPnlrowGetNLinearVars(nlrow)) );

   /* store all quadratic variables */
   SCIP_CALL( nlrowaggrStoreQuadraticVars(scip, *nlrowaggr, SCIPnlrowGetQuadVars(nlrow), SCIPnlrowGetNQuadVars(nlrow)) );
   assert((*nlrowaggr)->nquadvars == nquadvars);

   for( i = 0; i < nfound; ++i )
   {
      aggrnvars[i] = 0;
      aggrnterms[i] = 0;
   }

   /* count the number of variables in each e.c. aggregation */
   for( i = 0; i < nquadvars; ++i )
   {
      if( quadvar2aggr[i] >= 0)
         ++aggrnvars[ quadvar2aggr[i] ];
   }

   /* count the number of bilinear terms in each e.c. aggregation */
   for( i = 0; i < nquadelems; ++i )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_Real coef;
      int idx1;
      int idx2;

      quadelem = &SCIPnlrowGetQuadElems(nlrow)[i];
      idx1 = quadvar2aggr[quadelem->idx1];
      idx2 = quadvar2aggr[quadelem->idx2];
      coef = rhsaggr ? quadelem->coef : -quadelem->coef;

      /* variables has to belong to the same e.c. aggregation; bilinear / quadratic term has to be concave */
      if( idx1 >= 0 && idx2 >= 0 && idx1 == idx2 && (quadelem->idx1 != quadelem->idx2 || SCIPisNegative(scip, coef) ) )
         ++aggrnterms[idx1];
      else
         ++nremterms;
   }

   /* create all edge-concave aggregations (empty) and remaining terms */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlrowaggr)->ecaggr, nfound) );
   if( nremterms > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlrowaggr)->remtermcoefs, nremterms) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlrowaggr)->remtermvars1, nremterms) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlrowaggr)->remtermvars2, nremterms) );
      (*nlrowaggr)->remtermsize = nremterms;
   }
   (*nlrowaggr)->necaggr = nfound;

   for( i = 0; i < nfound; ++i )
   {
      SCIP_CALL( ecaggrCreateEmpty(scip, &(*nlrowaggr)->ecaggr[i], aggrnvars[i], aggrnterms[i]) );
   }

   /* add quadratic variables to the edge-concave aggregations */
   for( i = 0; i < nquadvars; ++i )
   {
      int idx;

      idx = quadvar2aggr[i];

      if( idx >= 0)
      {
         SCIPdebugMsg(scip, "add quadvar %d to aggr. %d\n", i, idx);
         SCIP_CALL( ecaggrAddQuadvar((*nlrowaggr)->ecaggr[idx], SCIPnlrowGetQuadVars(nlrow)[i]) );
      }
   }

   /* add the bilinear /quadratic terms to the edge-concave aggregations or in the remaining part */
   for( i = 0; i < nquadelems; ++i )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_Real coef;
      int idx1;
      int idx2;

      quadelem = &SCIPnlrowGetQuadElems(nlrow)[i];
      x = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1];
      y = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2];
      idx1 = quadvar2aggr[quadelem->idx1];
      idx2 = quadvar2aggr[quadelem->idx2];
      coef = rhsaggr ? quadelem->coef : -quadelem->coef;

      if( idx1 >= 0 && idx2 >= 0 && idx1 == idx2 && (quadelem->idx1 != quadelem->idx2 || SCIPisNegative(scip, coef) ) )
      {
         SCIP_CALL( ecaggrAddBilinTerm(scip, (*nlrowaggr)->ecaggr[idx1], x, y, coef) );
         SCIPdebugMsg(scip, "add term %e *%d*%d to aggr. %d\n", coef, quadelem->idx1, quadelem->idx2, idx1);
      }
      else
      {
         SCIP_CALL( nlrowaggrAddRemBilinTerm(scip, *nlrowaggr, x, y, coef) );
         SCIPdebugMsg(scip, "add term %e *%d*%d to the remaining part\n", coef, quadelem->idx1, quadelem->idx2);
      }
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &aggrnterms);
   SCIPfreeBufferArray(scip, &aggrnvars);

   return SCIP_OKAY;
}

/** frees a nonlinear row aggregation */
static
SCIP_RETCODE nlrowaggrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROWAGGR**      nlrowaggr           /**< pointer to free the nonlinear row aggregation */
   )
{
   int i;

   assert(scip != NULL);
   assert(nlrowaggr != NULL);
   assert(*nlrowaggr != NULL);
   (*nlrowaggr)->nlrow = NULL;
   assert((*nlrowaggr)->quadvars != NULL);
   assert((*nlrowaggr)->nquadvars > 0);
   assert((*nlrowaggr)->nremterms >= 0);

   /* free remaining part */
   SCIPfreeBlockMemoryArrayNull(scip, &(*nlrowaggr)->remtermcoefs, (*nlrowaggr)->remtermsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*nlrowaggr)->remtermvars1, (*nlrowaggr)->remtermsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*nlrowaggr)->remtermvars2, (*nlrowaggr)->remtermsize);

   /* free quadratic variables */
   SCIPfreeBlockMemoryArray(scip, &(*nlrowaggr)->quadvars, (*nlrowaggr)->nquadvars);
   SCIPfreeBlockMemoryArray(scip, &(*nlrowaggr)->quadvar2aggr, (*nlrowaggr)->nquadvars);
   (*nlrowaggr)->quadvars = NULL;
   (*nlrowaggr)->quadvar2aggr = NULL;
   (*nlrowaggr)->nquadvars = 0;

   /* free linear part */
   if( (*nlrowaggr)->nlinvars > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(*nlrowaggr)->linvars, (*nlrowaggr)->nlinvars);
      SCIPfreeBlockMemoryArray(scip, &(*nlrowaggr)->lincoefs, (*nlrowaggr)->nlinvars);
      (*nlrowaggr)->linvars = 0;
      (*nlrowaggr)->linvars = NULL;
      (*nlrowaggr)->lincoefs = NULL;
   }

   /* free edge-concave aggregations */
   for( i = 0; i < (*nlrowaggr)->necaggr; ++i )
   {
      SCIP_CALL( ecaggrFree(scip, &(*nlrowaggr)->ecaggr[i]) );
   }
   SCIPfreeBlockMemoryArray(scip, &(*nlrowaggr)->ecaggr, (*nlrowaggr)->necaggr);

   /* free nlrow aggregation */
   SCIPfreeBlockMemory(scip, nlrowaggr);

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** prints a nonlinear row aggregation */
static
void nlrowaggrPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROWAGGR*       nlrowaggr           /**< nonlinear row aggregation */
   )
{
   int i;

   SCIPdebugMsg(scip, " nlrowaggr rhs = %e\n", nlrowaggr->rhs);
   SCIPdebugMsg(scip, " #remaining terms = %d\n", nlrowaggr->nremterms);

   SCIPdebugMsg(scip, "remaining terms: ");
   for( i = 0; i < nlrowaggr->nremterms; ++i )
      SCIPdebugMsgPrint(scip, "%e %s * %s + ", nlrowaggr->remtermcoefs[i], SCIPvarGetName(nlrowaggr->remtermvars1[i]),
         SCIPvarGetName(nlrowaggr->remtermvars2[i]) );
   for( i = 0; i < nlrowaggr->nlinvars; ++i )
      SCIPdebugMsgPrint(scip, "%e %s + ", nlrowaggr->lincoefs[i], SCIPvarGetName(nlrowaggr->linvars[i]) );
   SCIPdebugMsgPrint(scip, "\n");

   for( i = 0; i < nlrowaggr->necaggr; ++i )
   {
      SCIPdebugMsg(scip, "print e.c. aggr %d\n", i);
      ecaggrPrint(scip, nlrowaggr->ecaggr[i]);
   }
   return;
}
#endif

/** creates separator data */
static
SCIP_RETCODE sepadataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< pointer to store separator data */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, sepadata) );
   BMSclearMemory(*sepadata);

   return SCIP_OKAY;
}

/** frees all nonlinear row aggregations */
static
SCIP_RETCODE sepadataFreeNlrows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< pointer to store separator data */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);

   /* free nonlinear row aggregations */
   if( sepadata->nlrowaggrs != NULL )
   {
      int i;

      for( i = sepadata->nnlrowaggrs - 1; i >= 0; --i )
      {
         SCIP_CALL( nlrowaggrFree(scip, &sepadata->nlrowaggrs[i]) );
      }

      SCIPfreeBlockMemoryArray(scip, &sepadata->nlrowaggrs, sepadata->nlrowaggrssize);

      sepadata->nlrowaggrs = NULL;
      sepadata->nnlrowaggrs = 0;
      sepadata->nlrowaggrssize = 0;
   }

   return SCIP_OKAY;
}

/** frees separator data */
static
SCIP_RETCODE sepadataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< pointer to store separator data */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(*sepadata != NULL);

   /* free nonlinear row aggregations */
   SCIP_CALL( sepadataFreeNlrows(scip, *sepadata) );

   /* free LP interface */
   if( (*sepadata)->lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&((*sepadata)->lpi)) );
      (*sepadata)->lpisize = 0;
   }

   SCIPfreeBlockMemory(scip, sepadata);

   return SCIP_OKAY;
}

/** adds a nonlinear row aggregation to the separator data */
static
SCIP_RETCODE sepadataAddNlrowaggr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROWAGGR*       nlrowaggr           /**< non-linear row aggregation */
   )
{
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(nlrowaggr != NULL);

   if( sepadata->nlrowaggrssize == 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->nlrowaggrs, 2) ); /*lint !e506*/
      sepadata->nlrowaggrssize = 2;
   }
   else if( sepadata->nlrowaggrssize < sepadata->nnlrowaggrs + 1 )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->nlrowaggrs, sepadata->nlrowaggrssize, 2 * sepadata->nlrowaggrssize) ); /*lint !e506 !e647*/
      sepadata->nlrowaggrssize *= 2;
      assert(sepadata->nlrowaggrssize >= sepadata->nnlrowaggrs + 1);
   }

   sepadata->nlrowaggrs[ sepadata->nnlrowaggrs ] = nlrowaggr;
   ++(sepadata->nnlrowaggrs);

   /* update maximum e.c. aggregation size */
   for( i = 0; i < nlrowaggr->necaggr; ++i )
      sepadata->maxecsize = MAX(sepadata->maxecsize, nlrowaggr->ecaggr[i]->nvars);

#ifdef SCIP_STATISTIC
   /* update statistics */
   if( nlrowaggr->rhsaggr )
      ++(sepadata->nrhsnlrowaggrs);
   else
      ++(sepadata->nlhsnlrowaggrs);
#endif

   return SCIP_OKAY;
}

/** returns min{val-lb,ub-val} / (ub-lb) */
static
SCIP_Real phi(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val,                /**< solution value */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   if( SCIPisFeasEQ(scip, lb, ub) )
      return 0.0;

   /* adjust */
   val = MAX(val, lb);
   val = MIN(val, ub);

   return MIN(ub - val, val - lb) / (ub - lb);
}

/** creates an MIP to search for cycles with an odd number of positive edges in the graph representation of a nonlinear
 *  row; the model uses directed binary arc flow variables; we introduce for all quadratic elements a forward and
 *  backward edge; if the term is quadratic (e.g., loop in the graph) we fix the corresponding variables to zero; this
 *  leads to an easy mapping of quadratic elements and the variables of the MIP
 */
static
SCIP_RETCODE createMIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< auxiliary SCIP to search aggregations */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_Bool             rhsaggr,            /**< consider nonlinear row aggregation for g(x) <= rhs (TRUE) or
                                              *   lhs <= g(x) (FALSE) */
   SCIP_VAR**            forwardarcs,        /**< array to store all forward arc variables */
   SCIP_VAR**            backwardarcs,       /**< array to store all backward arc variables */
   SCIP_Real*            nodeweights,        /**< weights for each node of the graph */
   int*                  nedges              /**< pointer to store the number of nonexcluded edges in the graph */
   )
{
   SCIP_VAR** oddcyclearcs;
   SCIP_CONS** flowcons;
   SCIP_CONS* cyclelengthcons;
   SCIP_CONS* oddcyclecons;
   char name[SCIP_MAXSTRLEN];
   int noddcyclearcs;
   int nnodes;
   int narcs;
   int i;

   assert(subscip != NULL);
   assert(forwardarcs != NULL);
   assert(backwardarcs != NULL);
   assert(nedges != NULL);
   assert(sepadata->minaggrsize <= sepadata->maxaggrsize);

   narcs = SCIPnlrowGetNQuadElems(nlrow);
   nnodes = SCIPnlrowGetNQuadVars(nlrow);
   *nedges = 0;

   assert(narcs > 0);
   assert(nnodes > 0);

   noddcyclearcs = 0;
   SCIP_CALL( SCIPallocBufferArray(subscip, &oddcyclearcs, 2*narcs) );

   /* create problem with default plug-ins */
   SCIP_CALL( SCIPcreateProbBasic(subscip, "E.C. aggregation MIP") );
   SCIP_CALL( SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* create forward and backward arc variables; loops are fixed to zero */
   for( i = 0; i < narcs; ++i )
   {
      SCIP_CONS* noparallelcons;
      SCIP_QUADELEM* quadelem;
      SCIP_Real edgeweight;
      SCIP_Real ub;

      quadelem  = &SCIPnlrowGetQuadElems(nlrow)[i];

      edgeweight = (quadelem->idx1 == quadelem->idx2) ? 0.0 : nodeweights[quadelem->idx1] + nodeweights[quadelem->idx2];
      SCIPdebugMsg(scip, "edge {%d,%d} = {%s,%s}  coeff=%e edgeweight=%e\n", quadelem->idx1, quadelem->idx2,
         SCIPvarGetName(SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1]),
         SCIPvarGetName(SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2]), quadelem->coef, edgeweight);

      ub = (quadelem->idx1 == quadelem->idx2) ? 0.0 : 1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x#%d#%d", quadelem->idx1, quadelem->idx2);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &forwardarcs[i], name, 0.0, ub, 0.01 + edgeweight, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(subscip, forwardarcs[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x#%d#%d", quadelem->idx2, quadelem->idx1);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &backwardarcs[i], name, 0.0, ub, 0.01 + edgeweight, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(subscip, backwardarcs[i]) );

      /* do not create redundant constraints for loops */
      if( quadelem->idx1 == quadelem->idx2 )
         continue;

      ++(*nedges);

      /* store all arcs which are important for the odd cycle property (no loops) */
      if( rhsaggr && SCIPisPositive(scip, quadelem->coef) )
      {
         oddcyclearcs[noddcyclearcs++] = forwardarcs[i];
         oddcyclearcs[noddcyclearcs++] = backwardarcs[i];
      }

      if( !rhsaggr && SCIPisNegative(scip, quadelem->coef) )
      {
         oddcyclearcs[noddcyclearcs++] = forwardarcs[i];
         oddcyclearcs[noddcyclearcs++] = backwardarcs[i];
      }

      /* add constraints to ensure no parallel edges  */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons_noparalleledges");
      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &noparallelcons, name, 0, NULL, NULL, 0.0, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, noparallelcons, forwardarcs[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, noparallelcons, backwardarcs[i], 1.0) );
      SCIP_CALL( SCIPaddCons(subscip, noparallelcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &noparallelcons) );
   }

   /* odd cycle property constraint */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons_oddcycle");
   SCIP_CALL( SCIPcreateConsBasicXor(subscip, &oddcyclecons, name, TRUE, noddcyclearcs, oddcyclearcs) );
   SCIP_CALL( SCIPaddCons(subscip, oddcyclecons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &oddcyclecons) );
   SCIPfreeBufferArray(subscip, &oddcyclearcs);

   /* cycle length constraint */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons_cyclelength");
   SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cyclelengthcons, name, 0, NULL, NULL,
         (SCIP_Real) sepadata->minaggrsize, (SCIP_Real) sepadata->maxaggrsize) );

   for( i = 0; i < narcs; ++i )
   {
      SCIP_CALL( SCIPaddCoefLinear(subscip, cyclelengthcons, forwardarcs[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, cyclelengthcons, backwardarcs[i], 1.0) );
   }

   SCIP_CALL( SCIPaddCons(subscip, cyclelengthcons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cyclelengthcons) );

   /* create flow conservation constraints */
   SCIP_CALL( SCIPallocBufferArray(subscip, &flowcons, nnodes) );

   for( i = 0; i < nnodes; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons_flowconservation#%d", i);
      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &flowcons[i], name, 0, NULL, NULL, 0.0, 0.0) );
   }

   for( i = 0; i < narcs; ++i )
   {
      int u;
      int v;

      u = SCIPnlrowGetQuadElems(nlrow)[i].idx1;
      assert(u >= 0 && u < SCIPnlrowGetNQuadVars(nlrow));
      v = SCIPnlrowGetQuadElems(nlrow)[i].idx2;
      assert(v >= 0 && v < SCIPnlrowGetNQuadVars(nlrow));

      SCIP_CALL( SCIPaddCoefLinear(subscip, flowcons[u], forwardarcs[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, flowcons[u], backwardarcs[i], -1.0) );

      SCIP_CALL( SCIPaddCoefLinear(subscip, flowcons[v], forwardarcs[i], -1.0) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, flowcons[v], backwardarcs[i], 1.0) );
   }

   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPaddCons(subscip, flowcons[i]) );
      SCIP_CALL( SCIPreleaseCons(subscip, &flowcons[i]) );
   }

   SCIPfreeBufferArray(subscip, &flowcons);

   return SCIP_OKAY;
}

/** fixed all arc variables (u,v) for which u or v is already in an edge-concave aggregation */
static
SCIP_RETCODE updateMIP(
   SCIP*                 subscip,            /**< auxiliary SCIP to search aggregations */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_VAR**            forwardarcs,        /**< forward arc variables */
   SCIP_VAR**            backwardarcs,       /**< backward arc variables */
   int*                  quadvar2aggr,       /**< mapping of quadvars to e.c. aggr. index (< 0: in no aggr.) */
   int*                  nedges              /**< pointer to store the number of nonexcluded edges */
   )
{
   int i;

   assert(subscip != NULL);
   assert(nlrow != NULL);
   assert(forwardarcs != NULL);
   assert(backwardarcs != NULL);
   assert(quadvar2aggr != NULL);
   assert(nedges != NULL);

   SCIP_CALL( SCIPfreeTransform(subscip) );

   /* recompute the number of edges */
   *nedges = 0;

   /* fix each arc to 0 if at least one of its nodes is contained in an e.c. aggregation */
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); ++i )
   {
      SCIP_QUADELEM* quadelem;

      quadelem = &SCIPnlrowGetQuadElems(nlrow)[i];

      if( quadvar2aggr[quadelem->idx1] != -1 || quadvar2aggr[quadelem->idx2] != -1 )
      {
         SCIP_CALL( SCIPchgVarUb(subscip, forwardarcs[i], 0.0) );
         SCIP_CALL( SCIPchgVarUb(subscip, backwardarcs[i], 0.0) );
      }
      else
         *nedges += (quadelem->idx1 != quadelem->idx2) ? 1 : 0;
   }

   return SCIP_OKAY;
}

/** stores the best edge-concave aggregation found by the MIP model */
static
SCIP_RETCODE storeAggrFromMIP(
   SCIP*                 subscip,            /**< auxiliary SCIP to search aggregations */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_VAR**            forwardarcs,        /**< forward arc variables */
   SCIP_VAR**            backwardarcs,       /**< backward arc variables */
   int*                  quadvar2aggr,       /**< mapping of quadvars to e.c. aggr. index (< 0: in no aggr.) */
   int                   nfoundsofar         /**< number of e.c. aggregation found so far */
   )
{
   SCIP_SOL* sol;
   int i;

   assert(subscip != NULL);
   assert(nlrow != NULL);
   assert(forwardarcs != NULL);
   assert(backwardarcs != NULL);
   assert(quadvar2aggr != NULL);
   assert(nfoundsofar >= 0);
   assert(SCIPgetStatus(subscip) < SCIP_STATUS_INFEASIBLE);
   assert(SCIPgetNSols(subscip) > 0);

   sol = SCIPgetBestSol(subscip);
   assert(sol != NULL);

   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); ++i )
   {
      SCIP_QUADELEM* quadelem;

      quadelem = &SCIPnlrowGetQuadElems(nlrow)[i];

      if( SCIPisGT(subscip, SCIPgetSolVal(subscip, sol, forwardarcs[i]), 0.5)
         || SCIPisGT(subscip, SCIPgetSolVal(subscip, sol, backwardarcs[i]), 0.5) )
      {
         assert(quadvar2aggr[quadelem->idx1] == -1 || quadvar2aggr[quadelem->idx1] == nfoundsofar);
         assert(quadvar2aggr[quadelem->idx2] == -1 || quadvar2aggr[quadelem->idx2] == nfoundsofar);

         quadvar2aggr[quadelem->idx1] = nfoundsofar;
         quadvar2aggr[quadelem->idx2] = nfoundsofar;
      }
   }

   return SCIP_OKAY;
}

/** searches for edge-concave aggregations with an MIP model based on binary flow variables */
static
SCIP_RETCODE searchEcAggrWithMIP(
   SCIP*                 subscip,            /**< SCIP data structure */
   SCIP_Real             timelimit,          /**< time limit to solve the MIP */
   int                   nedges,             /**< number of nonexcluded undirected edges */
   SCIP_Bool*            aggrleft,           /**< pointer to store if there might be a left aggregation */
   SCIP_Bool*            found               /**< pointer to store if we have found an aggregation */
   )
{
   assert(subscip != NULL);
   assert(aggrleft != NULL);
   assert(found != NULL);
   assert(nedges >= 0);

   *aggrleft = TRUE;
   *found = FALSE;

   if( SCIPisLE(subscip, timelimit, 0.0) )
      return SCIP_OKAY;

   /* set working limits */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/totalnodes", SUBSCIP_NODELIMIT) );

   /* set heuristics to aggressive */
   SCIP_CALL( SCIPsetHeuristics(subscip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );

   /* disable output to console in optimized mode, enable in SCIP's debug mode */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1) );
#else
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

   SCIP_CALL( SCIPsolve(subscip) );

   /* no more aggregation left if the MIP is infeasible */
   if( SCIPgetStatus(subscip) >= SCIP_STATUS_INFEASIBLE )
   {
      *found = FALSE;
      *aggrleft = FALSE;
      return SCIP_OKAY;
   }

   if( SCIPgetNSols(subscip) > 0 )
   {
      *found = TRUE;
      *aggrleft = TRUE;

#ifdef SCIP_DEBUG
      if( SCIPgetNSols(subscip) > 0 )
      {
         SCIP_CALL( SCIPprintSol(subscip, SCIPgetBestSol(subscip), NULL , FALSE) );
      }
#endif
   }

   return SCIP_OKAY;
}

/** creates a tclique graph from a given nonlinear row
 *
 *  SCIP's clique code can only handle integer node weights; all node weights are scaled by a factor of 100; since the
 *  clique code ignores nodes with weight of zero, we add an offset of 100 to each weight
 */
static
SCIP_RETCODE createTcliqueGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row  */
   TCLIQUE_GRAPH**       graph,              /**< TCLIQUE graph structure */
   SCIP_Real*            nodeweights         /**< weights for each quadratic variable (nodes in the graph) */
   )
{
   int i;

   assert(graph != NULL);
   assert(nlrow != NULL);

   /* create the tclique graph */
   if( !tcliqueCreate(graph) )
   {
      SCIPerrorMessage("could not create clique graph\n");
      return SCIP_ERROR;
   }

   /* add all nodes to the tclique graph */
   for( i = 0; i < SCIPnlrowGetNQuadVars(nlrow); ++i )
   {
      int nodeweight;

      /* note: clique code can only handle integer weights */
      nodeweight = 100 + (int)(100 * nodeweights[i]);
      /* SCIPdebugMsg(scip, "%d (%s): nodeweight %d \n", i, SCIPvarGetName(SCIPnlrowGetQuadVars(nlrow)[i]), nodeweight); */

      if( !tcliqueAddNode(*graph, i, nodeweight) )
      {
         SCIPerrorMessage("could not add node to clique graph\n");
         return SCIP_ERROR;
      }
   }

   /* add all edges */
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); ++i )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_VAR* bilinvar1;
      SCIP_VAR* bilinvar2;

      quadelem = &SCIPnlrowGetQuadElems(nlrow)[i];
      assert(quadelem != NULL);
      assert(!SCIPisZero(scip, quadelem->coef));

      bilinvar1 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1];
      bilinvar2 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2];
      assert(bilinvar1 != NULL);
      assert(bilinvar2 != NULL);

      /* do not add the edge {i,i} */
      if( bilinvar1 == bilinvar2 )
         continue;

      assert(quadelem->idx1 != quadelem->idx2);

#ifdef SCIP_DEBUG_DETAILED
      SCIPdebugMsg(scip, "   add edge (%d, %d) = (%s,%s) to tclique graph\n",
            SCIPvarGetIndex(bilinvar1), SCIPvarGetIndex(bilinvar2),
            SCIPvarGetName(bilinvar1), SCIPvarGetName(bilinvar2));
#endif

      if( !tcliqueAddEdge(*graph, quadelem->idx1, quadelem->idx2) )
      {
         SCIPerrorMessage("could not add edge to clique graph\n");
         return SCIP_ERROR;
      }
   }

   /* flush the clique graph */
   if( !tcliqueFlush(*graph) )
   {
      SCIPerrorMessage("could not flush the clique graph\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** searches for edge-concave aggregations by computing cliques in the graph representation of a given nonlinear row
 *  update graph, compute clique, store clique; after computing a clique we heuristically check if the clique contains
 *  at least one good cycle
 */
static
SCIP_RETCODE searchEcAggrWithCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,              /**< TCLIQUE graph structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   int*                  quadvar2aggr,       /**< mapping of quadvars to e.c. aggr. index (< 0: in no aggr.) */
   int                   nfoundsofar,        /**< number of e.c. aggregation found so far */
   SCIP_Bool             rhsaggr,            /**< consider nonlinear row aggregation for g(x) <= rhs (TRUE) or
                                              *   lhs <= g(x) (FALSE) */
   SCIP_Bool*            foundaggr,          /**< pointer to store if we have found an aggregation */
   SCIP_Bool*            foundclique         /**< pointer to store if we have found a clique */
   )
{
   SCIP_HASHMAP* cliquemap;
   TCLIQUE_STATUS status;
   int* maxcliquenodes;
   int* degrees;
   int nmaxcliquenodes;
   int maxcliqueweight;
   int noddcycleedges;
   int ntwodegrees;
   int aggrsize;
   int i;

   assert(graph != NULL);
   assert(nfoundsofar >= 0);
   assert(foundaggr != NULL);
   assert(foundclique != NULL);
   assert(SCIPnlrowGetNQuadVars(nlrow) == tcliqueGetNNodes(graph));

   cliquemap = NULL;
   *foundaggr = FALSE;
   *foundclique = FALSE;

   /* exclude all nodes which are already in an edge-concave aggregation (no flush is needed) */
   for( i = 0; i < SCIPnlrowGetNQuadVars(nlrow); ++i )
   {
      if( quadvar2aggr[i] != -1 )
      {
         SCIPdebugMsg(scip, "exclude node %d from clique graph\n", i);
         tcliqueChangeWeight(graph, i, 0);
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &maxcliquenodes, SCIPnlrowGetNQuadVars(nlrow)) );

   /* solve clique problem */
   tcliqueMaxClique(tcliqueGetNNodes, tcliqueGetWeights, tcliqueIsEdge, tcliqueSelectAdjnodes, graph, NULL, NULL,
      maxcliquenodes, &nmaxcliquenodes, &maxcliqueweight, CLIQUE_MAXFIRSTNODEWEIGHT, CLIQUE_MINWEIGHT,
      CLIQUE_MAXNTREENODES, CLIQUE_BACKTRACKFREQ, 0, -1, NULL, &status);

   if( status != TCLIQUE_OPTIMAL || nmaxcliquenodes < sepadata->minaggrsize )
      goto TERMINATE;

   *foundclique = TRUE;
   aggrsize = MIN(sepadata->maxaggrsize, nmaxcliquenodes);
   SCIP_CALL( SCIPhashmapCreate(&cliquemap, SCIPblkmem(scip), aggrsize) );

   for( i = 0; i < aggrsize; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(cliquemap, (void*) (size_t) maxcliquenodes[i], NULL) );
   }

   /* count the degree of good cycle edges for each node in the clique */
   SCIP_CALL( SCIPallocBufferArray(scip, &degrees, aggrsize) );
   BMSclearMemoryArray(degrees, aggrsize);
   ntwodegrees = 0;

   /* count the number of positive or negative edges (depending on <= rhs or >= lhs) */
   noddcycleedges = 0;
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); ++i )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_Bool isoddcycleedge;

      quadelem = &SCIPnlrowGetQuadElems(nlrow)[i];
      isoddcycleedge = (rhsaggr && SCIPisPositive(scip, quadelem->coef))
         || (!rhsaggr && SCIPisNegative(scip, quadelem->coef));

      if( isoddcycleedge
         && SCIPhashmapExists(cliquemap, (void*) (size_t) quadelem->idx1)
         && SCIPhashmapExists(cliquemap, (void*) (size_t) quadelem->idx2) )
      {
         ++noddcycleedges;
         ++degrees[quadelem->idx1];
         ++degrees[quadelem->idx2];
      }
   }

   /* count the number of nodes with exactly two incident odd cycle edges */
   for( i = 0; i < aggrsize; ++i )
      if( degrees[i] == 2 )
         ++ntwodegrees;

   /* check cases for which we are sure that there are no good cycles in the clique */
   if( noddcycleedges == 0 || (aggrsize == 3 && noddcycleedges == 2) || (aggrsize == 4 && ntwodegrees == 4) )
      *foundaggr = FALSE;
   else
      *foundaggr = TRUE;

   /* add the found clique as an edge-concave aggregation or exclude the nodes from the remaining search */
   for( i = 0; i < aggrsize; ++i )
   {
      quadvar2aggr[ maxcliquenodes[i] ] = *foundaggr ? nfoundsofar : -2;
      SCIPdebugMsg(scip, "%s %d\n", *foundaggr ? "aggregate node: " : "exclude node: ", maxcliquenodes[i]);
   }

   SCIPfreeBufferArray(scip, &degrees);

TERMINATE:
   if( cliquemap != NULL )
      SCIPhashmapFree(&cliquemap);
   SCIPfreeBufferArray(scip, &maxcliquenodes);

   return SCIP_OKAY;
}

/** helper function for searchEcAggr() */
static
SCIP_RETCODE doSeachEcAggr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SOL*             sol,                /**< current solution (might be NULL) */
   SCIP_Bool             rhsaggr,            /**< consider nonlinear row aggregation for g(x) <= rhs (TRUE) or g(x) >= lhs (FALSE) */
   int*                  quadvar2aggr,       /**< array to store for each quadratic variable in which edge-concave
                                              *   aggregation it is stored (< 0: in no aggregation); size has to be at
                                              *   least SCIPnlrowGetNQuadVars(nlrow) */
   int*                  nfound              /**< pointer to store the number of found e.c. aggregations */
   )
{
   TCLIQUE_GRAPH* graph = NULL;
   SCIP_VAR** forwardarcs;
   SCIP_VAR** backwardarcs;
   SCIP_VAR** quadvars;
   SCIP_Real* nodeweights;
   SCIP_Real timelimit;
   int nunsucces = 0;
   int nedges = 0;
   int nquadelems;
   int nquadvars;
   int i;

   assert(subscip != NULL);
   assert(quadvar2aggr != NULL);
   assert(nfound != NULL);

   quadvars = SCIPnlrowGetQuadVars(nlrow);
   nquadvars = SCIPnlrowGetNQuadVars(nlrow);
   nquadelems = SCIPnlrowGetNQuadElems(nlrow);

   *nfound = 0;

   /* arrays to store all arc variables of the MIP model; note that we introduce variables even for loops in the graph
    * to have an easy mapping from the edges of the graph to the quadratic elements
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeweights, nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &forwardarcs, nquadelems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &backwardarcs, nquadelems) );

   /* inititialize mapping from quadvars to e.c. aggregation index (-1: quadvar is in no aggregation); compute node
    * weights
    */
   for( i = 0; i < SCIPnlrowGetNQuadVars(nlrow); ++i )
   {
      SCIP_VAR* var = quadvars[i];
      quadvar2aggr[i] = -1;
      nodeweights[i] = phi(scip, SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var) );
      SCIPdebugMsg(scip, "%s = %e (%e in [%e, %e])\n", SCIPvarGetName(var), nodeweights[i], SCIPgetSolVal(scip, sol, var),
         SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
   }

   SCIP_CALL( createMIP(scip, subscip, sepadata, nlrow, rhsaggr, forwardarcs, backwardarcs, nodeweights, &nedges) );
   assert(nedges >= 0);
   SCIPdebugMsg(scip, "nedges (without loops) = %d\n", nedges);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* main loop to search for edge-concave aggregations */
   while( !SCIPisStopped(scip) )
   {
      SCIP_Bool aggrleft;
      SCIP_Bool found;

      SCIPdebugMsg(scip, "#remaining edges = %d\n", nedges);

      /* not enough edges left */
      if( nedges < sepadata->minaggrsize )
         break;

      /* check whether there is enough time left; update the remaining time */
      if( !SCIPisInfinity(scip, timelimit) )
      {
         timelimit -= SCIPgetSolvingTime(scip);
         if( timelimit <= 0.0 )
         {
            SCIPdebugMsg(scip, "skip aggregation search since no time left\n");
            goto TERMINATE;
         }
      }

      /* 1.a - search for edge-concave aggregation with the help of the MIP model */
      SCIP_CALL( searchEcAggrWithMIP(subscip, timelimit, nedges, &aggrleft, &found) );

      /* 1.b - there are no more edge-concave aggregations left */
      if( !aggrleft )
      {
         SCIPdebugMsg(scip, "no more aggregation left\n");
         break;
      }

      if( found )
      {
         SCIP_CALL( storeAggrFromMIP(subscip, nlrow, forwardarcs, backwardarcs, quadvar2aggr, *nfound) );
         ++(*nfound);
         nunsucces = 0;
      }
      /* try to find an edge-concave aggregation by computing cliques */
      else
      {
         SCIP_Bool foundaggr;
         SCIP_Bool foundclique;

         ++nunsucces;

         /* create graph if necessary */
         if( graph == NULL )
         {
            SCIP_CALL( createTcliqueGraph(scip, nlrow, &graph, nodeweights) );
            assert(graph != NULL);
         }

         /* 2.a - search and store a single edge-concave aggregation by computing a clique with a good cycle */
         SCIP_CALL_FINALLY( searchEcAggrWithCliques(scip, graph, sepadata, nlrow, quadvar2aggr, *nfound, rhsaggr,
               &foundaggr, &foundclique), tcliqueFree(&graph) );

         if( foundaggr )
         {
            assert(foundclique);
            ++(*nfound);
            nunsucces = 0;
         }
         else
            ++nunsucces;

         /* 2.b - no clique of at least minaggrsize size found */
         if( !foundclique )
         {
            assert(!foundaggr);
            SCIPdebugMsg(scip, "did not find a clique to exclude -> leave aggregation search\n");
            break;
         }
      }

      /* leave the algorithm if we did not find something for maxstallrounds many iterations */
      if( nunsucces >= sepadata->maxstallrounds && *nfound == 0 )
      {
         SCIPdebugMsg(scip, "did not find an e.c. aggregation for %d iterations\n", nunsucces);
         break;
      }

      /* exclude all edges used in the last aggregation and nodes found in the clique solution */
      SCIP_CALL_FINALLY( updateMIP(subscip, nlrow, forwardarcs, backwardarcs, quadvar2aggr, &nedges), tcliqueFree(&graph) );
   }

TERMINATE:

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "aggregations found:\n");
   for( i = 0; i < nquadvars; ++i )
   {
      SCIPdebugMsg(scip, " %d in %d\n", i, quadvar2aggr[i]);
   }
#endif

   /* free clique graph */
   if( graph != NULL )
      tcliqueFree(&graph);

   /* free sub-SCIP */
   for( i = 0; i < nquadelems; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &forwardarcs[i]) );
      SCIP_CALL( SCIPreleaseVar(subscip, &backwardarcs[i]) );
   }

   SCIPfreeBufferArray(scip, &backwardarcs);
   SCIPfreeBufferArray(scip, &forwardarcs);
   SCIPfreeBufferArray(scip, &nodeweights);

   return SCIP_OKAY;
}

/** computes a partitioning into edge-concave aggregations for a given (quadratic) nonlinear row; each aggregation has
 *  to contain a cycle with an odd number of positive weighted edges (good cycles) in the corresponding graph representation
 *
 *  For this we use the following algorithm:
 *
 *  -# use a MIP model based on binary flow variables to compute good cycles and store the implied subgraphs as an e.c. aggr.
 *    -# if we find a good cycle, store the implied subgraph, delete it from the graph representation and go to 1)
 *    -# if the MIP model is infeasible (there are no good cycles), STOP
 *  -# we compute a large clique C if the MIP model fails (because of working limits, etc)
 *    -# if we find a good cycle in C, store the implied subgraph of C, delete it from the graph representation and go to 1)
 *    -# if C is not large enough, STOP
 */
static
SCIP_RETCODE searchEcAggr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SOL*             sol,                /**< current solution (might be NULL) */
   SCIP_Bool             rhsaggr,            /**< consider nonlinear row aggregation for g(x) <= rhs (TRUE) or g(x) >= lhs (FALSE) */
   int*                  quadvar2aggr,       /**< array to store for each quadratic variable in which edge-concave
                                              *   aggregation it is stored (< 0: in no aggregation); size has to be at
                                              *   least SCIPnlrowGetNQuadVars(nlrow) */
   int*                  nfound              /**< pointer to store the number of found e.c. aggregations */
   )
{
   SCIP* subscip;
   SCIP_RETCODE retcode;

   /* create and set up a sub-SCIP */
   SCIP_CALL_FINALLY( SCIPcreate(&subscip), (void)SCIPfree(&subscip) );

   retcode = doSeachEcAggr(scip, subscip, sepadata, nlrow, sol, rhsaggr, quadvar2aggr, nfound);

   SCIP_CALL( SCIPfree(&subscip) );
   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** returns whether a given nonlinear row can be used to compute edge-concave aggregations for which their convex
 *  envelope could dominate the termwise bilinear relaxation; this is the case if there exists at least one cycle with
 *  an odd number of positive edges in the corresponding graph representation of the nonlinear row
 */
static
SCIP_RETCODE isCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW*           nlrow,              /**< nonlinear row representation of a nonlinear constraint */
   SCIP_Bool*            rhscandidate,       /**< pointer to store if we should compute edge-concave aggregations for
                                              *   the <= rhs case */
   SCIP_Bool*            lhscandidate        /**< pointer to store if we should compute edge-concave aggregations for
                                              *   the >= lhs case */
   )
{
   int* degrees;
   int ninterestingnodes;
   int nposedges;
   int nnegedges;
   int i;

   assert(rhscandidate != NULL);
   assert(lhscandidate != NULL);

   *rhscandidate = TRUE;
   *lhscandidate = TRUE;

   /* skip if the nlrow is not in the NLP, there are other nonlinearities, or too few quadratic variables */
   if( !SCIPnlrowIsInNLP(nlrow) || SCIPnlrowGetExprtree(nlrow) != NULL
      || SCIPnlrowGetNQuadVars(nlrow) < sepadata->minaggrsize )
   {
      *rhscandidate = FALSE;
      *lhscandidate = FALSE;
      return SCIP_OKAY;
   }

   /* check for infinite rhs or lhs */
   if( SCIPisInfinity(scip, REALABS(SCIPnlrowGetRhs(nlrow))) )
      *rhscandidate = FALSE;
   if( SCIPisInfinity(scip, REALABS(SCIPnlrowGetLhs(nlrow))) )
      *lhscandidate = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &degrees, SCIPnlrowGetNQuadVars(nlrow)) );
   BMSclearMemoryArray(degrees, SCIPnlrowGetNQuadVars(nlrow));

   ninterestingnodes = 0;
   nposedges = 0;
   nnegedges = 0;

   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); ++i )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_VAR* x;
      SCIP_VAR* y;

      quadelem = &SCIPnlrowGetQuadElems(nlrow)[i];
      x = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1];
      y = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2];

      /* do not consider loops or global fixed variables */
      if( quadelem->idx1 != quadelem->idx2
         && !SCIPisEQ(scip, SCIPvarGetLbGlobal(x), SCIPvarGetUbGlobal(x))
         && !SCIPisEQ(scip, SCIPvarGetLbGlobal(y), SCIPvarGetUbGlobal(y)) )
      {
         ++degrees[quadelem->idx1];
         ++degrees[quadelem->idx2];

         /* count the number of nodes with a degree of at least 2 */
         if( degrees[quadelem->idx1] == 2 )
            ++ninterestingnodes;
         if( degrees[quadelem->idx2] == 2 )
            ++ninterestingnodes;

         nposedges += SCIPisPositive(scip, quadelem->coef) ? 1 : 0;
         nnegedges += SCIPisNegative(scip, quadelem->coef) ? 1 : 0;
      }
   }

   SCIPfreeBufferArray(scip, &degrees);

   SCIPdebugMsg(scip, "nlrow contains: %d edges\n", nposedges + nnegedges);

   /* too many edges, too few edges, or to few nodes with degree at least 2 in the graph  */
   if( nposedges + nnegedges > sepadata->maxbilinterms || nposedges + nnegedges < sepadata->minaggrsize
      || ninterestingnodes < sepadata->minaggrsize )
   {
      *rhscandidate = FALSE;
      *lhscandidate = FALSE;
      return SCIP_OKAY;
   }

   /* check if there are enough positive/negative edges; for a 3-clique there has to be an odd number of those edges */
   if( nposedges == 0 || (nposedges + nnegedges == 3 && (nposedges % 2) == 0) )
      *rhscandidate = FALSE;
   if( nnegedges == 0 || (nposedges + nnegedges == 3 && (nnegedges % 2) == 0) )
      *lhscandidate = FALSE;

   return SCIP_OKAY;
}

/** finds and stores edge-concave aggregations for a given nonlinear row */
static
SCIP_RETCODE findAndStoreEcAggregations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SOL*             sol                 /**< current solution (might be NULL) */
   )
{
   int* quadvar2aggr;
   SCIP_Bool rhscandidate;
   SCIP_Bool lhscandidate;

   assert(scip != NULL);
   assert(nlrow != NULL);
   assert(sepadata != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &quadvar2aggr, SCIPnlrowGetNQuadVars(nlrow)) ); /*lint !e705*/

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "search for edge-concave aggregation for the nonlinear row: \n");
   SCIP_CALL( SCIPnlrowPrint(nlrow, SCIPgetMessagehdlr(scip), NULL) );
#endif

   /* check obvious conditions for existing cycles with an odd number of positive/negative edges */
   SCIP_CALL( isCandidate(scip, sepadata, nlrow, &rhscandidate, &lhscandidate) );
   SCIPdebugMsg(scip, "rhs candidate = %u lhs candidate = %u\n", rhscandidate, lhscandidate);

   /* search for edge-concave aggregations (consider <= rhs) */
   if( rhscandidate )
   {
      SCIP_NLROWAGGR* nlrowaggr;
      int nfound;

      assert(!SCIPisInfinity(scip, REALABS(SCIPnlrowGetRhs(nlrow))));

      SCIPdebugMsg(scip, "consider <= rhs\n");
      SCIP_CALL( searchEcAggr(scip, sepadata, nlrow, sol, TRUE, quadvar2aggr, &nfound) );

      if( nfound > 0 )
      {
         SCIP_CALL( nlrowaggrCreate(scip, nlrow, &nlrowaggr, quadvar2aggr, nfound, TRUE) );
         assert(nlrow != NULL);
         SCIPdebug(nlrowaggrPrint(scip, nlrowaggr));
         SCIP_CALL( sepadataAddNlrowaggr(scip, sepadata, nlrowaggr) );
      }
   }

   /* search for edge-concave aggregations (consider <= lhs) */
   if( lhscandidate )
   {
      SCIP_NLROWAGGR* nlrowaggr;
      int nfound;

      assert(!SCIPisInfinity(scip, REALABS(SCIPnlrowGetLhs(nlrow))));

      SCIPdebugMsg(scip, "consider >= lhs\n");
      SCIP_CALL( searchEcAggr(scip, sepadata, nlrow, sol, FALSE, quadvar2aggr, &nfound) );

      if( nfound > 0 )
      {
         SCIP_CALL( nlrowaggrCreate(scip, nlrow, &nlrowaggr, quadvar2aggr, nfound, FALSE) );
         assert(nlrow != NULL);
         SCIPdebug(nlrowaggrPrint(scip, nlrowaggr));
         SCIP_CALL( sepadataAddNlrowaggr(scip, sepadata, nlrowaggr) );
      }
   }

   SCIPfreeBufferArray(scip, &quadvar2aggr);
   return SCIP_OKAY;
}

/*
 *  methods to compute edge-concave cuts
 */

#ifdef SCIP_DEBUG
/** prints a given facet (candidate) */
static
void printFacet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables contained in the edge-concave aggregation */
   int                   nvars,              /**< number of variables contained in the edge-concave aggregation */
   SCIP_Real*            facet,              /**< current facet candidate */
   SCIP_Real             facetval            /**< facet evaluated at the current solution */
   )
{
   int i;

   SCIPdebugMsg(scip, "print facet (val=%e): ", facetval);
   for( i = 0; i < nvars; ++i )
      SCIPdebugMsgPrint(scip, "%e %s + ", facet[i], SCIPvarGetName(vars[i]));
   SCIPdebugMsgPrint(scip, "%e\n", facet[nvars]);
}
#endif

/** checks if a facet is really an underestimate for all corners of the domain [l,u]; because of numerics it can happen
 *  that a facet violates a corner of the domain; to make the facet valid we subtract the maximum violation from the
 *  constant part of the facet; its a good excersise to write a comment describing the gray code...
 */
static
SCIP_Bool checkRikun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ECAGGR*          ecaggr,             /**< edge-concave aggregation data */
   SCIP_Real*            fvals,              /**< array containing all corner values of the aggregation */
   SCIP_Real*            facet               /**< current facet candidate (of dimension ecaggr->nvars + 1) */
   )
{
   SCIP_Real maxviolation;
   SCIP_Real val;
   unsigned int i;
   unsigned int ncorner;
   unsigned int prev;

   assert(scip != NULL);
   assert(ecaggr != NULL);
   assert(fvals != NULL);
   assert(facet != NULL);

   ncorner = (unsigned int) poweroftwo[ecaggr->nvars];
   maxviolation = 0.0;

   /* check for the origin */
   val = facet[ecaggr->nvars];
   for( i = 0; i < (unsigned int) ecaggr->nvars; ++i )
      val += facet[i] * SCIPvarGetLbLocal(ecaggr->vars[i]);

   /* update  maximum violation */
   maxviolation = MAX(val - fvals[0], maxviolation);
   assert(SCIPisFeasEQ(scip, maxviolation, 0.0));

   prev = 0;
   for( i = 1; i < ncorner; ++i )
   {
      unsigned int gray;
      unsigned int diff;
      unsigned int pos;

      gray = i ^ (i >> 1);
      diff = gray ^ prev;

      /* compute position of unique 1 of diff */
      pos = 0;
      while( (diff >>= 1) != 0 )
         ++pos;

      if( gray > prev )
         val += facet[pos] * (SCIPvarGetUbLocal(ecaggr->vars[pos]) - SCIPvarGetLbLocal(ecaggr->vars[pos]));
      else
         val -= facet[pos] * (SCIPvarGetUbLocal(ecaggr->vars[pos]) - SCIPvarGetLbLocal(ecaggr->vars[pos]));


      /* update  maximum violation */
      maxviolation = MAX(val - fvals[gray], maxviolation);
      assert(SCIPisFeasEQ(scip, maxviolation, 0.0));

      prev = gray;
   }

   SCIPdebugMsg(scip, "maximum violation of facet: %2.8e\n", maxviolation);

   /* there seem to be numerical problems if the violation is too large; in this case we reject the facet */
   if( maxviolation > ADJUSTFACETTOL )
      return FALSE;

   /* adjust constant part of the facet */
   facet[ecaggr->nvars] -= maxviolation;

   return TRUE;
}

/** set up LP interface to solve LPs to compute the facet of the convex envelope */
static
SCIP_RETCODE createLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separation data */
   )
{
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* val;
   int* beg;
   int* ind;
   int nnonz;
   int ncols;
   int nrows;
   int i;
   int k;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->nnlrowaggrs > 0);

   /* LP interface has been already created with enough rows/columns*/
   if( sepadata->lpi != NULL && sepadata->lpisize >= sepadata->maxecsize )
      return SCIP_OKAY;

   /* size of lpi is too small; reconstruct lpi */
   if( sepadata->lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&sepadata->lpi) );
      sepadata->lpi = NULL;
   }

   assert(sepadata->lpi == NULL);
   SCIP_CALL( SCIPlpiCreate(&(sepadata->lpi), SCIPgetMessagehdlr(scip), "e.c. LP", SCIP_OBJSEN_MINIMIZE) );
   sepadata->lpisize = sepadata->maxecsize;

   nrows = sepadata->maxecsize + 1;
   ncols = poweroftwo[nrows - 1];
   nnonz = (ncols * (nrows + 1)) / 2;
   k = 0;

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, nnonz) );

   /* calculate nonzero entries in the LP; set obj, lb, and ub to zero */
   for( i = 0; i < ncols; ++i )
   {
      int row;
      int a;

      obj[i] = 0.0;
      lb[i] = 0.0;
      ub[i] = 0.0;

      SCIPdebugMsg(scip, "col %i starts at position %d\n", i, k);
      beg[i] = k;
      row = 0;
      a = 1;

      /* iterate through the bit representation of i */
      while( a <= i )
      {
         if( (a & i) != 0 )
         {
            val[k] = 1.0;
            ind[k] = row;

            SCIPdebugMsg(scip, " val[%d][%d] = 1 (position  %d)\n", row, i, k);

            ++k;
         }

         a <<= 1; /*lint !e701*/
         ++row;
         assert(poweroftwo[row] == a);
      }

      /* put 1 as a coefficient for sum_{i} \lambda_i = 1 row (last row) */
      val[k] = 1.0;
      ind[k] = nrows - 1;
      ++k;
      SCIPdebugMsg(scip, " val[%d][%d] = 1 (position  %d)\n", nrows - 1, i, k);
   }
   assert(k == nnonz);

   /*
    * add all columns to the LP interface
    * CPLEX needs the row to exist before adding columns, so we create the rows with dummy sides
    * note that the assert is not needed once somebody fixes the LPI
    */
   assert(nrows <= ncols);
   SCIP_CALL( SCIPlpiAddRows(sepadata->lpi, nrows, obj, obj, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(sepadata->lpi, ncols, obj, lb, ub, NULL, nnonz, beg, ind, val) );

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** evaluates an edge-concave aggregation at a corner of the domain [l,u] */
static
SCIP_Real evalCorner(
   SCIP_ECAGGR*          ecaggr,             /**< edge-concave aggregation data */
   int                   k                   /**< k-th corner */
   )
{
   SCIP_Real val;
   int i;

   assert(ecaggr != NULL);
   assert(k >= 0 && k < poweroftwo[ecaggr->nvars]);

   val = 0.0;

   for( i = 0; i < ecaggr->nterms; ++i )
   {
      SCIP_Real coef;
      SCIP_Real bound1;
      SCIP_Real bound2;
      int idx1;
      int idx2;

      idx1 = ecaggr->termvars1[i];
      idx2 = ecaggr->termvars2[i];
      coef = ecaggr->termcoefs[i];
      assert(idx1 >= 0 && idx1 < ecaggr->nvars);
      assert(idx2 >= 0 && idx2 < ecaggr->nvars);

      bound1 = ((poweroftwo[idx1]) & k) == 0 ? SCIPvarGetLbLocal(ecaggr->vars[idx1]) : SCIPvarGetUbLocal(ecaggr->vars[idx1]);
      bound2 = ((poweroftwo[idx2]) & k) == 0 ? SCIPvarGetLbLocal(ecaggr->vars[idx2]) : SCIPvarGetUbLocal(ecaggr->vars[idx2]);

      val += coef * bound1 * bound2;
   }

   return val;
}

/** returns (val - lb) / (ub - lb) for a in [lb, ub] */
static
SCIP_Real transformValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub,                 /**< upper bound */
   SCIP_Real             val                 /**< value in [lb,ub] */
   )
{
   assert(scip != NULL);
   assert(!SCIPisInfinity(scip, -lb));
   assert(!SCIPisInfinity(scip, ub));
   assert(!SCIPisInfinity(scip, REALABS(val)));
   assert(!SCIPisFeasEQ(scip, ub - lb, 0.0)); /* this would mean that a variable has been fixed */

   /* adjust val */
   val = MIN(val, ub);
   val = MAX(val, lb);

   val = (val - lb) / (ub - lb);
   assert(val >= 0.0 && val <= 1.0);

   return val;
}

/** computes a facet of the convex envelope of an edge concave aggregation
 *
 *  The algorithm solves the following LP:
 *  \f{eqnarray}{
 *      min   & \sum_i \lambda_i f(v_i)\\
 *      s.t.  & \sum_i \lambda_i v_i = x\\
 *            & \sum_i \lambda_i     = 1\\
 *            & \lambda_i \geq 0
 *  \f}
 *  where f is an edge concave function, \f$x\f$ in \f$[l,u]\f$ is a solution of the current relaxation, and \f$v_i\f$ are the vertices
 *  of \f$[l,u]\f$; the method transforms the problem to the domain \f$[0,1]^n\f$, computes a facet, and transforms this facet to the
 *  original space; the dual solution of the LP above are the coefficients of the facet
 *
 *  The complete algorithm works as follows:
 *
 *  -# compute f(v_i) for each corner v_i of [l,u]
 *  -# set up the described LP for the transformed space
 *  -# solve the LP and store the resulting facet for the transformed space
 *  -# transform the facet to original space
 *  -# adjust and check facet with the algorithm of Rikun et al.
 */
static
SCIP_RETCODE SCIPcomputeConvexEnvelopeFacet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separation data */
   SCIP_SOL*             sol,                /**< solution (might be NULL) */
   SCIP_ECAGGR*          ecaggr,             /**< edge-concave aggregation data */
   SCIP_Real*            facet,              /**< array to store the coefficients of the resulting facet; size has to be at least (ecaggr->nvars + 1) */
   SCIP_Real*            facetval,           /**< pointer to store the value of the facet evaluated at the current solution */
   SCIP_Bool*            success             /**< pointer to store if we have found a facet */
   )
{
   SCIP_Real* fvals;
   SCIP_Real* side;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real perturbation;
   int* inds;
   int ncorner;
   int ncols;
   int nrows;
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(ecaggr != NULL);
   assert(facet != NULL);
   assert(facetval != NULL);
   assert(success != NULL);
   assert(ecaggr->nvars <= sepadata->maxecsize);

   *facetval = -SCIPinfinity(scip);
   *success = FALSE;

   /* create LP if this has not been done yet */
   SCIP_CALL( createLP(scip, sepadata) );

   assert(sepadata->lpi != NULL);
   assert(sepadata->lpisize >= ecaggr->nvars);

   SCIP_CALL( SCIPlpiGetNCols(sepadata->lpi, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(sepadata->lpi, &nrows) );
   ncorner = poweroftwo[ecaggr->nvars];

   assert(ncorner <= ncols);
   assert(ecaggr->nvars + 1 <= nrows);
   assert(nrows <= ncols);

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &fvals, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &side, ncols) );

   /*
    *  1. compute f(v_i) for each corner v_i of [l,u]
    *  2. set up the described LP for the transformed space
    */
   for( i = 0; i < ncols; ++i )
   {
      fvals[i] = i < ncorner ? evalCorner(ecaggr, i) : 0.0;
      inds[i] = i;

      /* update bounds; fix variables to zero which are currently not in the LP */
      lb[i] = 0.0;
      ub[i] = i < ncorner ? 1.0 : 0.0;
      SCIPdebugMsg(scip, "bounds of LP col %d = [%e, %e]; obj = %e\n", i, lb[i], ub[i], fvals[i]);
   }

   /* update lhs and rhs */
   perturbation = 0.001;
   for( i = 0; i < nrows; ++i )
   {
      /* note that the last row corresponds to sum_{j} \lambda_j = 1 */
      if( i < ecaggr->nvars )
      {
         SCIP_VAR* x;

         x = ecaggr->vars[i];
         assert(x != NULL);

         side[i] = transformValue(scip, SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), SCIPgetSolVal(scip, sol, x));

         /* perturb point to enforce an LP solution with ecaggr->nvars + 1 nonzero */
         side[i] += side[i] > perturbation ? -perturbation : perturbation;
         perturbation /= 1.2;
      }
      else
      {
         side[i] = (i == nrows - 1) ? 1.0 : 0.0;
      }

      SCIPdebugMsg(scip, "LP row %d in [%e, %e]\n", i, side[i], side[i]);
   }

   /* update LP */
   SCIP_CALL( SCIPlpiChgObj(sepadata->lpi, ncols, inds, fvals) );
   SCIP_CALL( SCIPlpiChgBounds(sepadata->lpi, ncols, inds, lb, ub) );
   SCIP_CALL( SCIPlpiChgSides(sepadata->lpi, nrows, inds, side, side) );

   /* free memory used to build the LP */
   SCIPfreeBufferArray(scip, &side);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &inds);

   /*
    *  3. solve the LP and store the resulting facet for the transformed space
    */
   if( USEDUALSIMPLEX ) /*lint !e774 !e506*/
   {
      SCIP_CALL( SCIPlpiSolveDual(sepadata->lpi) );
   }
   else
   {
      SCIP_CALL( SCIPlpiSolvePrimal(sepadata->lpi) );
   }

   /* the dual solution corresponds to the coefficients of the facet in the transformed problem; note that it might be
    * the case that the dual solution has more components than the facet array
    */
   if( ecaggr->nvars + 1 == ncols )
   {
      SCIP_CALL( SCIPlpiGetSol(sepadata->lpi, NULL, NULL, facet, NULL, NULL) );
   }
   else
   {
      SCIP_Real* dualsol;

      SCIP_CALL( SCIPallocBufferArray(scip, &dualsol, nrows) );

      /* get the dual solution */
      SCIP_CALL( SCIPlpiGetSol(sepadata->lpi, NULL, NULL, dualsol, NULL, NULL) );

      for( i = 0; i < ecaggr->nvars; ++i )
         facet[i] = dualsol[i];

      /* constant part of the facet is the last component of the dual solution */
      facet[ecaggr->nvars] = dualsol[nrows - 1];

      SCIPfreeBufferArray(scip, &dualsol);
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "facet for the transformed problem: ");
   for( i = 0; i < ecaggr->nvars; ++i )
   {
      SCIPdebugMsgPrint(scip, "%3.4e * %s + ", facet[i], SCIPvarGetName(ecaggr->vars[i]));
   }
   SCIPdebugMsgPrint(scip, "%3.4e\n", facet[ecaggr->nvars]);
#endif

   /*
    *  4. transform the facet to original space
    *  we now have the linear underestimator L(x) = beta^T x + beta_0, which needs to be transform to the original space
    *  the underestimator in the original space, G(x) = alpha^T x + alpha_0, is given by G(x) = L(T(x)), where T(.) is
    *  the transformation applied in step 2; therefore,
    *  alpha_i = beta_i/(ub_i - lb_i)
    *  alpha_0 = beta_0 - sum_i lb_i * beta_i/(ub_i - lb_i)
    */

   SCIPdebugMsg(scip, "facet in orig. space: ");
   *facetval = 0.0;

   for( i = 0; i < ecaggr->nvars; ++i )
   {
      SCIP_Real varlb;
      SCIP_Real varub;

      varlb = SCIPvarGetLbLocal(ecaggr->vars[i]);
      varub = SCIPvarGetUbLocal(ecaggr->vars[i]);
      assert(!SCIPisEQ(scip, varlb, varub));

      /* substract (\beta_i * lb_i) / (ub_i - lb_i) from current alpha_0 */
      facet[ecaggr->nvars] -= (facet[i] * varlb) / (varub - varlb);

      /* set \alpha_i := \beta_i / (ub_i - lb_i) */
      facet[i] = facet[i] / (varub - varlb);
      *facetval += facet[i] * SCIPgetSolVal(scip, sol, ecaggr->vars[i]);

      SCIPdebugMsgPrint(scip, "%3.4e * %s + ", facet[i], SCIPvarGetName(ecaggr->vars[i]));
   }

   /* add constant part to the facet value */
   *facetval += facet[ecaggr->nvars];
   SCIPdebugMsgPrint(scip, "%3.4e\n", facet[ecaggr->nvars]);

   /*
    *  5. adjust and check facet with the algorithm of Rikun et al.
    */

   if( checkRikun(scip, ecaggr, fvals, facet) )
   {
      SCIPdebugMsg(scip, "facet pass the check of Rikun et al.\n");
      *success = TRUE;
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &fvals);

   return SCIP_OKAY;
}

/*
 * miscellaneous methods
 */

/** method to add a facet of the convex envelope of an edge-concave aggregation to a given cut */
static
SCIP_RETCODE addFacetToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< current solution (might be NULL) */
   SCIP_ROW*             cut,                /**< current cut (modifiable) */
   SCIP_Real*            facet,              /**< coefficient of the facet (dimension nvars + 1) */
   SCIP_VAR**            vars,               /**< variables of the facet */
   int                   nvars,              /**< number of variables in the facet */
   SCIP_Real*            cutconstant,        /**< pointer to update the constant part of the facet */
   SCIP_Real*            cutactivity,        /**< pointer to update the activity of the cut */
   SCIP_Bool*            success             /**< pointer to store if everything went fine */
   )
{
   int i;

   assert(cut != NULL);
   assert(facet != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(cutconstant != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);

   *success = TRUE;

   for( i = 0; i < nvars; ++i )
   {
      if( SCIPisInfinity(scip, REALABS(facet[i])) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      if( !SCIPisZero(scip, facet[i]) )
      {
         /* add only a constant if the variable has been fixed */
         if( SCIPvarGetLbLocal(vars[i]) == SCIPvarGetUbLocal(vars[i]) ) /*lint !e777*/
         {
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[i]), SCIPgetSolVal(scip, sol, vars[i])));
            *cutconstant += facet[i] * SCIPgetSolVal(scip, sol, vars[i]);
            *cutactivity += facet[i] * SCIPgetSolVal(scip, sol, vars[i]);
         }
         else
         {
            *cutactivity += facet[i] * SCIPgetSolVal(scip, sol, vars[i]);
            SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[i], facet[i]) );
         }
      }
   }

   /* add constant part of the facet */
   *cutconstant += facet[nvars];
   *cutactivity += facet[nvars];

   return SCIP_OKAY;
}

/** method to add an linear term to a given cut */
static
SCIP_RETCODE addLinearTermToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< current solution (might be NULL) */
   SCIP_ROW*             cut,                /**< current cut (modifiable) */
   SCIP_VAR*             x,                  /**< linear variable */
   SCIP_Real             coeff,              /**< coefficient */
   SCIP_Real*            cutconstant,        /**< pointer to update the constant part of the facet */
   SCIP_Real*            cutactivity,        /**< pointer to update the activity of the cut */
   SCIP_Bool*            success             /**< pointer to store if everything went fine */
   )
{
   SCIP_Real activity;

   assert(cut != NULL);
   assert(x != NULL);
   assert(!SCIPisZero(scip, coeff));
   assert(!SCIPisInfinity(scip, coeff));
   assert(cutconstant != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);

   *success = TRUE;
   activity = SCIPgetSolVal(scip, sol, x) * coeff;

   /* do not add a term if the activity is -infinity */
   if( SCIPisInfinity(scip, -1.0 * REALABS(activity)) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* add activity to the constant part if the variable has been fixed */
   if( SCIPvarGetLbLocal(x) == SCIPvarGetUbLocal(x) ) /*lint !e777*/
   {
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(x), SCIPgetSolVal(scip, sol, x)));
      *cutconstant += activity;
      SCIPdebugMsg(scip, "add to cut: %e\n", activity);
   }
   else
   {
      SCIP_CALL( SCIPaddVarToRow(scip, cut, x, coeff) );
      SCIPdebugMsg(scip, "add to cut: %e * %s\n", coeff, SCIPvarGetName(x));
   }

   *cutactivity += activity;

   return SCIP_OKAY;
}

/** method to add an underestimate of a bilinear term to a given cut */
static
SCIP_RETCODE addBilinearTermToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< current solution (might be NULL) */
   SCIP_ROW*             cut,                /**< current cut (modifiable) */
   SCIP_VAR*             x,                  /**< first bilinear variable */
   SCIP_VAR*             y,                  /**< seconds bilinear variable */
   SCIP_Real             coeff,              /**< coefficient */
   SCIP_Real*            cutconstant,        /**< pointer to update the constant part of the facet */
   SCIP_Real*            cutactivity,        /**< pointer to update the activity of the cut */
   SCIP_Bool*            success             /**< pointer to store if everything went fine */
   )
{
   SCIP_Real activity;

   assert(cut != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(!SCIPisZero(scip, coeff));
   assert(cutconstant != NULL);
   assert(cutactivity != NULL);
   assert(success != NULL);

   *success = TRUE;
   activity = coeff * SCIPgetSolVal(scip, sol, x) * SCIPgetSolVal(scip, sol, y);

   if( SCIPisInfinity(scip, REALABS(coeff)) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* do not add a term if the activity is -infinity */
   if( SCIPisInfinity(scip, -1.0 * REALABS(activity)) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* quadratic case */
   if( x == y )
   {
      SCIP_Real refpoint;
      SCIP_Real lincoef;
      SCIP_Real linconst;

      lincoef = 0.0;
      linconst = 0.0;
      refpoint = SCIPgetSolVal(scip, sol, x);

      /* adjust the reference point */
      refpoint = SCIPisLT(scip, refpoint, SCIPvarGetLbLocal(x)) ? SCIPvarGetLbLocal(x) : refpoint;
      refpoint = SCIPisGT(scip, refpoint, SCIPvarGetUbLocal(x)) ? SCIPvarGetUbLocal(x) : refpoint;
      assert(SCIPisLE(scip, refpoint, SCIPvarGetUbLocal(x)) && SCIPisGE(scip, refpoint, SCIPvarGetLbLocal(x)));

      if( SCIPisPositive(scip, coeff) )
         SCIPaddSquareLinearization(scip, coeff, refpoint, SCIPvarIsIntegral(x), &lincoef, &linconst, success);
      else
         SCIPaddSquareSecant(scip, coeff, SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), refpoint, &lincoef, &linconst, success);

      *cutactivity += lincoef * refpoint + linconst;
      *cutconstant += linconst;

      /* add underestimate to cut */
      SCIP_CALL( SCIPaddVarToRow(scip, cut, x, lincoef) );

      SCIPdebugMsg(scip, "add to cut: %e * %s + %e\n", lincoef, SCIPvarGetName(x), linconst);
   }
   /* bilinear case */
   else
   {
      SCIP_Real refpointx;
      SCIP_Real refpointy;
      SCIP_Real lincoefx;
      SCIP_Real lincoefy;
      SCIP_Real linconst;

      lincoefx = 0.0;
      lincoefy = 0.0;
      linconst = 0.0;
      refpointx = SCIPgetSolVal(scip, sol, x);
      refpointy = SCIPgetSolVal(scip, sol, y);

      /* adjust the reference points */
      refpointx = SCIPisLT(scip, refpointx, SCIPvarGetLbLocal(x)) ? SCIPvarGetLbLocal(x) : refpointx;
      refpointx = SCIPisGT(scip, refpointx, SCIPvarGetUbLocal(x)) ? SCIPvarGetUbLocal(x) : refpointx;
      refpointy = SCIPisLT(scip, refpointy, SCIPvarGetLbLocal(y)) ? SCIPvarGetLbLocal(y) : refpointy;
      refpointy = SCIPisGT(scip, refpointy, SCIPvarGetUbLocal(y)) ? SCIPvarGetUbLocal(y) : refpointy;
      assert(SCIPisLE(scip, refpointx, SCIPvarGetUbLocal(x)) && SCIPisGE(scip, refpointx, SCIPvarGetLbLocal(x)));
      assert(SCIPisLE(scip, refpointy, SCIPvarGetUbLocal(y)) && SCIPisGE(scip, refpointy, SCIPvarGetLbLocal(y)));

      SCIPaddBilinMcCormick(scip, coeff, SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), refpointx, SCIPvarGetLbLocal(y),
         SCIPvarGetUbLocal(y), refpointy, FALSE, &lincoefx, &lincoefy, &linconst, success);

      *cutactivity += lincoefx * refpointx + lincoefy * refpointy + linconst;
      *cutconstant += linconst;

      /* add underestimate to cut */
      SCIP_CALL( SCIPaddVarToRow(scip, cut, x, lincoefx) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, y, lincoefy) );

      SCIPdebugMsg(scip, "add to cut: %e * %s + %e * %s + %e\n", lincoefx, SCIPvarGetName(x), lincoefy,
         SCIPvarGetName(y), linconst);
   }

   return SCIP_OKAY;
}

/** method to compute and and a cut for a nonlinear row aggregation and a given solution; we compute for each edge
 *  concave aggregation one facet; the remaining bilinear terms will be underestimated with McCormick, secants or
 *  linearizations; constant and linear terms will be added to the cut directly
 */
static
SCIP_RETCODE computeCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROWAGGR*       nlrowaggr,          /**< nonlinear row aggregation */
   SCIP_SOL*             sol,                /**< current solution (might be NULL) */
   SCIP_Bool*            separated,          /**< pointer to store if we could separate the current solution */
   SCIP_Bool*            cutoff              /**< pointer to store if the current node gets cut off */
   )
{
   SCIP_ROW* cut;
   SCIP_Real* bestfacet;
   SCIP_Real bestfacetval;
   SCIP_Real cutconstant;
   SCIP_Real cutactivity;
   int bestfacetsize;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Bool found;
   SCIP_Bool islocalcut;
   int i;

   assert(separated != NULL);
   assert(cutoff != NULL);
   assert(nlrowaggr->necaggr > 0);
   assert(nlrowaggr->nlrow != NULL);
   assert(SCIPnlrowIsInNLP(nlrowaggr->nlrow));

   *separated = FALSE;
   *cutoff = FALSE;
   islocalcut = SCIPgetDepth(scip) != 0;

   /* create the cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "ec");
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), SCIPinfinity(scip), islocalcut, FALSE,
         sepadata->dynamiccuts) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

   /* track rhs and activity of the cut */
   cutconstant = nlrowaggr->constant;
   cutactivity = 0.0;

   /* allocate necessary memory */
   bestfacetsize = sepadata->maxaggrsize + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &bestfacet, bestfacetsize) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPnlrowPrint(nlrowaggr->nlrow, SCIPgetMessagehdlr(scip), NULL) );

   SCIPdebugMsg(scip, "current solution:\n");
   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      SCIP_VAR* var = SCIPgetVars(scip)[i];
      SCIPdebugMsg(scip, "  %s = [%e, %e] solval = %e\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
         SCIPvarGetUbLocal(var), SCIPgetSolVal(scip, sol, var));
   }
#endif

   /* compute a facet for each edge-concave aggregation */
   for( i = 0; i < nlrowaggr->necaggr; ++i )
   {
      SCIP_ECAGGR* ecaggr;
      SCIP_Bool success;

      ecaggr = nlrowaggr->ecaggr[i];
      assert(ecaggr != NULL);

      /* compute a facet of the convex envelope */
      SCIP_CALL( SCIPcomputeConvexEnvelopeFacet(scip, sepadata, sol, ecaggr, bestfacet, &bestfacetval, &found) );

      SCIPdebugMsg(scip, "found facet for edge-concave aggregation %d/%d ? %s\n", i, nlrowaggr->necaggr,
         found ? "yes" : "no");

#ifdef SCIP_DEBUG
      if( found )
         printFacet(scip, ecaggr->vars, ecaggr->nvars, bestfacet, bestfacetval);
#endif

      /* do not add any cut because we did not found a facet for at least one edge-concave aggregation */
      if( !found ) /*lint !e774*/
         goto TERMINATE;

      /* add facet to the cut and update the rhs and activity of the cut */
      SCIP_CALL( addFacetToCut(scip, sol, cut, bestfacet, ecaggr->vars, ecaggr->nvars, &cutconstant, &cutactivity,
            &success) );

      if( !success )
         goto TERMINATE;
   }

   /* compute an underestimate for each bilinear term which is not in any edge-concave aggregation */
   for( i = 0; i < nlrowaggr->nremterms; ++i )
   {
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_Bool success;

      x = nlrowaggr->remtermvars1[i];
      y = nlrowaggr->remtermvars2[i];
      assert(x != NULL);
      assert(y != NULL);

      SCIP_CALL( addBilinearTermToCut(scip, sol, cut, x, y, nlrowaggr->remtermcoefs[i], &cutconstant, &cutactivity,
            &success) );

      if( !success )
         goto TERMINATE;
   }

   /* add all linear terms to the cut */
   for( i = 0; i < nlrowaggr->nlinvars; ++i )
   {
      SCIP_VAR* x;
      SCIP_Real coef;
      SCIP_Bool success;

      x = nlrowaggr->linvars[i];
      assert(x != NULL);

      coef = nlrowaggr->lincoefs[i];

      SCIP_CALL( addLinearTermToCut(scip, sol, cut, x, coef, &cutconstant, &cutactivity, &success) );

      if( !success )
         goto TERMINATE;
   }

   SCIPdebugMsg(scip, "cut activity = %e rhs(nlrow) = %e\n", cutactivity, nlrowaggr->rhs);

   /* set rhs of the cut (substract the constant part of the cut) */
   SCIP_CALL( SCIPchgRowRhs(scip, cut, nlrowaggr->rhs - cutconstant) );
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

   /* check activity of the row; this assert can fail because of numerics */
   /* assert(SCIPisFeasEQ(scip, cutactivity - cutconstant, SCIPgetRowSolActivity(scip, cut, sol)) ); */

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

   SCIPdebugMsg(scip, "EC cut <%s>: act=%f eff=%f rank=%d range=%e\n",
      SCIProwGetName(cut), SCIPgetRowSolActivity(scip, cut, sol), SCIPgetCutEfficacy(scip, sol, cut),
      SCIProwGetRank(cut), SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut) );

   /* try to add the cut has a finite rhs, is efficacious, and does not exceed the maximum cut range */
   if( !SCIPisInfinity(scip, nlrowaggr->rhs - cutconstant) && SCIPisCutEfficacious(scip, sol, cut)
      && SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut) < sepadata->cutmaxrange )
   {
      /* add the cut if it is separating the given solution by at least minviolation */
      if( SCIPisGE(scip, cutactivity - nlrowaggr->rhs, sepadata->minviolation) )
      {
         SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
         *separated = TRUE;
         SCIPdebugMsg(scip, "added separating cut\n");
      }

      if( !(*cutoff) && !islocalcut )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         SCIPdebugMsg(scip, "added cut to cut pool\n");
      }
   }

TERMINATE:
   /* free allocated memory */
   SCIPfreeBufferArray(scip, &bestfacet);

   /* release the row */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** returns whether it is possible to compute a cut for a given nonlinear row aggregation */
static
SCIP_Bool isPossibleToComputeCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< current solution (might be NULL) */
   SCIP_NLROWAGGR*       nlrowaggr           /**< nonlinear row aggregation */
   )
{
   int i;

   assert(scip != NULL);
   assert(nlrowaggr != NULL);

   if( !SCIPnlrowIsInNLP(nlrowaggr->nlrow) )
   {
      SCIPdebugMsg(scip, "nlrow is not in NLP anymore\n");
      return FALSE;
   }

   for( i = 0; i < nlrowaggr->nquadvars; ++i )
   {
      SCIP_VAR* var = nlrowaggr->quadvars[i];
      assert(var != NULL);

      /* check whether the variable has infinite bounds */
      if( SCIPisInfinity(scip, REALABS(SCIPvarGetLbLocal(var))) || SCIPisInfinity(scip, REALABS(SCIPvarGetUbLocal(var)))
            || SCIPisInfinity(scip, REALABS(SCIPgetSolVal(scip, sol, var))) )
      {
         SCIPdebugMsg(scip, "nlrow aggregation contains unbounded variables\n");
         return FALSE;
      }

      /* check whether the variable has been fixed and is in one edge-concave aggregation */
      if( nlrowaggr->quadvar2aggr[i] >= 0 && SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMsg(scip, "nlrow aggregation contains fixed variables in an e.c. aggregation\n");
         return FALSE;
      }
   }

   return TRUE;
}

/** searches and tries to add edge-concave cuts */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< current solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   int nmaxcuts;
   int ncuts;
   int i;

   assert(*result == SCIP_DIDNOTRUN);

   SCIPdebugMsg(scip, "separate cuts...\n");

   /* skip if there are no nonlinear row aggregations */
   if( sepadata->nnlrowaggrs == 0 )
   {
      SCIPdebugMsg(scip, "no aggregations exists -> skip call\n");
      return SCIP_OKAY;
   }

   /* get the maximal number of cuts allowed in a separation round */
   nmaxcuts = SCIPgetDepth(scip) == 0 ? sepadata->maxsepacutsroot : sepadata->maxsepacuts;
   ncuts = 0;

   /* try to compute cuts for each nonlinear row independently */
   for( i = 0; i < sepadata->nnlrowaggrs && ncuts < nmaxcuts && !SCIPisStopped(scip); ++i )
   {
      SCIP_NLROWAGGR* nlrowaggr;
      SCIP_Bool separated;
      SCIP_Bool cutoff;

      nlrowaggr = sepadata->nlrowaggrs[i];
      assert(nlrowaggr != NULL);

      /* skip nonlinear aggregations for which it is obviously not possible to compute a cut */
      if( !isPossibleToComputeCut(scip, sol, nlrowaggr) )
         return SCIP_OKAY;

      *result = (*result == SCIP_DIDNOTRUN) ? SCIP_DIDNOTFIND : *result;

      SCIPdebugMsg(scip, "try to compute a cut for nonlinear row aggregation %d\n", i);

      /* compute and add cut */
      SCIP_CALL( computeCut(scip, sepa, sepadata, nlrowaggr, sol, &separated, &cutoff) );
      SCIPdebugMsg(scip, "found a cut: %s cutoff: %s\n", separated ? "yes" : "no", cutoff ? "yes" : "no");

      /* stop if the current node gets cut off */
      if( cutoff )
      {
         assert(separated);
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* do not compute more cuts if we already separated the given solution */
      if( separated )
      {
         assert(!cutoff);
         *result = SCIP_SEPARATED;
         ++ncuts;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyEccuts)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaEccuts(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeEccuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( sepadataFree(scip, &sepadata) );
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolEccuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* print statistics */
#ifdef SCIP_STATISTIC
   SCIPstatisticMessage("rhs-AGGR %d\n", sepadata->nrhsnlrowaggrs);
   SCIPstatisticMessage("lhs-AGGR %d\n", sepadata->nlhsnlrowaggrs);
   SCIPstatisticMessage("aggr. search time = %f\n", sepadata->aggrsearchtime);
#endif

   /* free nonlinear row aggregations */
   SCIP_CALL( sepadataFreeNlrows(scip, sepadata) );

   /* mark that we should search again for nonlinear row aggregations */
   sepadata->searchedforaggr = FALSE;

   SCIPdebugMsg(scip, "exitsol\n");

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpEccuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   int depth;
   int ncalls;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( !allowlocal )
      return SCIP_OKAY;

   /* check min- and maximal aggregation size */
   if( sepadata->maxaggrsize < sepadata->minaggrsize )
      return SCIP_PARAMETERWRONGVAL;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* skip if the LP is not constructed yet */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "Skip since NLP is not constructed yet.\n");
      return SCIP_OKAY;
   }

   depth = SCIPgetDepth(scip);

   /* only call separator up to a maximum depth */
   if ( sepadata->maxdepth >= 0 && depth > sepadata->maxdepth )
      return SCIP_OKAY;

   /* only call separator a given number of times at each node */
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   if ( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* search for nonlinear row aggregations */
   if( !sepadata->searchedforaggr )
   {
      int i;

      SCIPstatistic( sepadata->aggrsearchtime -= SCIPgetTotalTime(scip) );

      SCIPdebugMsg(scip, "search for nonlinear row aggregations\n");
      for( i = 0; i < SCIPgetNNLPNlRows(scip) && !SCIPisStopped(scip); ++i )
      {
         SCIP_NLROW* nlrow = SCIPgetNLPNlRows(scip)[i];
         SCIP_CALL( findAndStoreEcAggregations(scip, sepadata, nlrow, NULL) );
      }
      sepadata->searchedforaggr = TRUE;

      SCIPstatistic( sepadata->aggrsearchtime += SCIPgetTotalTime(scip) );
   }

   /* search for edge-concave cuts */
   SCIP_CALL( separateCuts(scip, sepa, sepadata, NULL, result) );

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the edge concave separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaEccuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create eccuts separator data */
   SCIP_CALL( sepadataCreate(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpEccuts, NULL, sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyEccuts) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeEccuts) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolEccuts) );

   /* add eccuts separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxrounds",
         "maximal number of eccuts separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxroundsroot",
         "maximal number of eccuts separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxdepth",
         "maximal depth at which the separator is applied (-1: unlimited)",
         &sepadata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacuts",
         "maximal number of edge-concave cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacutsroot",
         "maximal number of edge-concave cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/cutmaxrange",
         "maximal coef. range of a cut (max coef. divided by min coef.) in order to be added to LP relaxation",
         &sepadata->cutmaxrange, FALSE, DEFAULT_CUTMAXRANGE, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/minviolation",
         "minimal violation of an edge-concave cut to be separated",
         &sepadata->minviolation, FALSE, DEFAULT_MINVIOLATION, 0.0, 0.5, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/minaggrsize",
         "search for edge-concave aggregations of at least this size",
         &sepadata->minaggrsize, TRUE, DEFAULT_MINAGGRSIZE, 3, 5, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxaggrsize",
         "search for edge-concave aggregations of at most this size",
         &sepadata->maxaggrsize, TRUE, DEFAULT_MAXAGGRSIZE, 3, 5, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxbilinterms",
         "maximum number of bilinear terms allowed to be in a quadratic constraint",
         &sepadata->maxbilinterms, TRUE, DEFAULT_MAXBILINTERMS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxstallrounds",
         "maximum number of unsuccessful rounds in the edge-concave aggregation search",
         &sepadata->maxstallrounds, TRUE, DEFAULT_MAXSTALLROUNDS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
