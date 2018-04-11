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

/**@file   cons_knapsack.c
 * @brief  Constraint handler for knapsack constraints of the form  \f$a^T x \le b\f$, x binary and \f$a \ge 0\f$.
 * @author Tobias Achterberg
 * @author Xin Liu
 * @author Kati Wolter
 * @author Michael Winkler
 * @author Tobias Fischer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <ctype.h>

#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/pub_misc.h"

#ifdef WITH_CARDINALITY_UPGRADE
#include "scip/cons_cardinality.h"
#endif

/* constraint handler properties */
#define CONSHDLR_NAME          "knapsack"
#define CONSHDLR_DESC          "knapsack constraint of the form  a^T x <= b, x binary and a >= 0"
#define CONSHDLR_SEPAPRIORITY   +600000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -600000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -600000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_ALWAYS
#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "knapsack"
#define EVENTHDLR_DESC         "bound change event handler for knapsack constraints"
#define EVENTTYPE_KNAPSACK SCIP_EVENTTYPE_LBCHANGED \
                         | SCIP_EVENTTYPE_UBTIGHTENED \
                         | SCIP_EVENTTYPE_VARFIXED \
                         | SCIP_EVENTTYPE_VARDELETED \
                         | SCIP_EVENTTYPE_IMPLADDED /**< variable events that should be caught by the event handler */

#define LINCONSUPGD_PRIORITY    +100000 /**< priority of the constraint handler for upgrading of linear constraints */

#define MAX_USECLIQUES_SIZE        1000 /**< maximal number of items in knapsack where clique information is used */
#define MAX_ZEROITEMS_SIZE        10000 /**< maximal number of items to store in the zero list in preprocessing */

#define KNAPSACKRELAX_MAXDELTA        0.1 /**< maximal allowed rounding distance for scaling in knapsack relaxation */
#define KNAPSACKRELAX_MAXDNOM      1000LL /**< maximal allowed denominator in knapsack rational relaxation */
#define KNAPSACKRELAX_MAXSCALE     1000.0 /**< maximal allowed scaling factor in knapsack rational relaxation */

#define DEFAULT_SEPACARDFREQ          1 /**< multiplier on separation frequency, how often knapsack cuts are separated */
#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of cuts separated per separation round in the root node */
#define DEFAULT_MAXCARDBOUNDDIST    0.0 /**< maximal relative distance from current node's dual bound to primal bound compared
                                         *   to best node's dual bound for separating knapsack cuts */
#define DEFAULT_DISAGGREGATION     TRUE /**< should disaggregation of knapsack constraints be allowed in preprocessing? */
#define DEFAULT_SIMPLIFYINEQUALITIES TRUE/**< should presolving try to simplify knapsacks */
#define DEFAULT_NEGATEDCLIQUE      TRUE /**< should negated clique information be used in solving process */

#define MAXABSVBCOEF               1e+5 /**< maximal absolute coefficient in variable bounds used for knapsack relaxation */
#define USESUPADDLIFT             FALSE /**< should lifted minimal cover inequalities using superadditive up-lifting be separated in addition */

#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define HASHSIZE_KNAPSACKCONS       500 /**< minimal size of hash table in linear constraint tables */

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presolving comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise 
                                         *   comparison round */
#define DEFAULT_DUALPRESOLVING     TRUE /**< should dual presolving steps be performed? */
#define DEFAULT_DETECTCUTOFFBOUND  TRUE /**< should presolving try to detect constraints parallel to the objective
                                         *   function defining an upper bound and prevent these constraints from
                                         *   entering the LP */
#define DEFAULT_DETECTLOWERBOUND TRUE   /**< should presolving try to detect constraints parallel to the objective
                                         *   function defining a lower bound and prevent these constraints from
                                         *   entering the LP */
#define DEFAULT_CLIQUEEXTRACTFACTOR 0.5 /**< lower clique size limit for greedy clique extraction algorithm (relative to largest clique) */
#define MAXCOVERSIZEITERLEWI       1000 /**< maximal size for which LEWI are iteratively separated by reducing the feasible set */

#define DEFAULT_USEGUBS           FALSE /**< should GUB information be used for separation? */
#define GUBCONSGROWVALUE              6 /**< memory growing value for GUB constraint array */
#define GUBSPLITGNC1GUBS          FALSE /**< should GNC1 GUB conss without F vars be split into GOC1 and GR GUB conss? */
#define DEFAULT_CLQPARTUPDATEFAC   1.5  /**< factor on the growth of global cliques to decide when to update a previous
                                         *   (negated) clique partition (used only if updatecliquepartitions is set to TRUE) */
#define DEFAULT_UPDATECLIQUEPARTITIONS FALSE /**< should clique partition information be updated when old partition seems outdated? */
#define MAXNCLIQUEVARSCOMP 1000000      /**< limit on number of pairwise comparisons in clique partitioning algorithm */
#ifdef WITH_CARDINALITY_UPGRADE
#define DEFAULT_UPGDCARDINALITY   FALSE /**< if TRUE then try to update knapsack constraints to cardinality constraints */
#endif

/* @todo maybe use event SCIP_EVENTTYPE_VARUNLOCKED to decide for another dual-presolving run on a constraint */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int*                  ints1;              /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   int*                  ints2;              /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Longint*         longints1;          /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Longint*         longints2;          /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools1;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools2;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools3;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Bool*            bools4;             /**< cleared memory array, all entries are set to zero in initpre, if you use this
                                              *   you have to clear it at the end, exists only in presolving stage */
   SCIP_Real*            reals1;             /**< cleared memory array, all entries are set to zero in consinit, if you use this
                                              *   you have to clear it at the end */
   int                   ints1size;          /**< size of ints1 array */
   int                   ints2size;          /**< size of ints2 array */
   int                   longints1size;      /**< size of longints1 array */
   int                   longints2size;      /**< size of longints2 array */
   int                   bools1size;         /**< size of bools1 array */
   int                   bools2size;         /**< size of bools2 array */
   int                   bools3size;         /**< size of bools3 array */
   int                   bools4size;         /**< size of bools4 array */
   int                   reals1size;         /**< size of reals1 array */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_Real             maxcardbounddist;   /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for separating knapsack cuts */
   int                   sepacardfreq;       /**< multiplier on separation frequency, how often knapsack cuts are separated */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in the root node */
   SCIP_Bool             disaggregation;     /**< should disaggregation of knapsack constraints be allowed in preprocessing? */
   SCIP_Bool             simplifyinequalities;/**< should presolving try to cancel down or delete coefficients in inequalities */
   SCIP_Bool             negatedclique;      /**< should negated clique information be used in solving process */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             dualpresolving;     /**< should dual presolving steps be performed? */
   SCIP_Bool             usegubs;            /**< should GUB information be used for separation? */
   SCIP_Bool             detectcutoffbound;  /**< should presolving try to detect constraints parallel to the objective
                                              *   function defining an upper bound and prevent these constraints from
                                              *   entering the LP */
   SCIP_Bool             detectlowerbound;   /**< should presolving try to detect constraints parallel to the objective
                                              *   function defining a lower bound and prevent these constraints from
                                              *   entering the LP */
   SCIP_Bool             updatecliquepartitions; /**< should clique partition information be updated when old partition seems outdated? */
   SCIP_Real             cliqueextractfactor;/**< lower clique size limit for greedy clique extraction algorithm (relative to largest clique) */
   SCIP_Real             clqpartupdatefac;   /**< factor on the growth of global cliques to decide when to update a previous
                                              *   (negated) clique partition (used only if updatecliquepartitions is set to TRUE) */
#ifdef WITH_CARDINALITY_UPGRADE
   SCIP_Bool             upgdcardinality;    /**< if TRUE then try to update knapsack constraints to cardinality constraints */
   SCIP_Bool             upgradedcard;       /**< whether we have already upgraded knapsack constraints to cardinality constraints */
#endif
};


/** constraint data for knapsack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in knapsack constraint */
   SCIP_Longint*         weights;            /**< weights of variables in knapsack constraint */
   SCIP_EVENTDATA**      eventdata;          /**< event data for bound change events of the variables */
   int*                  cliquepartition;    /**< clique indices of the clique partition */
   int*                  negcliquepartition; /**< clique indices of the negated clique partition */
   SCIP_ROW*             row;                /**< corresponding LP row */
   int                   nvars;              /**< number of variables in knapsack constraint */
   int                   varssize;           /**< size of vars, weights, and eventdata arrays */
   int                   ncliques;           /**< number of cliques in the clique partition */
   int                   nnegcliques;        /**< number of cliques in the negated clique partition */
   int                   ncliqueslastnegpart;/**< number of global cliques the last time a negated clique partition was computed */
   int                   ncliqueslastpart;   /**< number of global cliques the last time a clique partition was computed */
   SCIP_Longint          capacity;           /**< capacity of knapsack */
   SCIP_Longint          weightsum;          /**< sum of all weights */
   SCIP_Longint          onesweightsum;      /**< sum of weights of variables fixed to one */
   unsigned int          presolvedtiming:5;  /**< max level in which the knapsack constraint is already presolved */
   unsigned int          sorted:1;           /**< are the knapsack items sorted by weight? */
   unsigned int          cliquepartitioned:1;/**< is the clique partition valid? */
   unsigned int          negcliquepartitioned:1;/**< is the negated clique partition valid? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the knapsack already added to clique table? */
   unsigned int          varsdeleted:1;      /**< were variables deleted after last cleanup? */
   unsigned int          existmultaggr:1;    /**< does this constraint contain multi-aggregations */
};

/** event data for bound changes events */
struct SCIP_EventData
{
   SCIP_CONS*            cons;               /**< knapsack constraint to process the bound change for */
   SCIP_Longint          weight;             /**< weight of variable */
   int                   filterpos;          /**< position of event in variable's event filter */
};


/** data structure to combine two sorting key values */
struct sortkeypair
{
   SCIP_Real             key1;               /**< first sort key value */
   SCIP_Real             key2;               /**< second sort key value */
};
typedef struct sortkeypair SORTKEYPAIR;

/** status of GUB constraint */
enum GUBVarstatus
{
   GUBVARSTATUS_UNINITIAL        = -1,       /** unintitialized variable status */
   GUBVARSTATUS_CAPACITYEXCEEDED =  0,       /** variable with weight exceeding the knapsack capacity */
   GUBVARSTATUS_BELONGSTOSET_R   =  1,       /** variable in noncovervars R */
   GUBVARSTATUS_BELONGSTOSET_F   =  2,       /** variable in noncovervars F */
   GUBVARSTATUS_BELONGSTOSET_C2  =  3,       /** variable in covervars C2 */
   GUBVARSTATUS_BELONGSTOSET_C1  =  4        /** variable in covervars C1 */
};
typedef enum GUBVarstatus GUBVARSTATUS;

/** status of variable in GUB constraint */
enum GUBConsstatus
{
   GUBCONSSTATUS_UNINITIAL         = -1,     /** unintitialized GUB constraint status */
   GUBCONSSTATUS_BELONGSTOSET_GR   =  0,     /** all GUB variables are in noncovervars R */
   GUBCONSSTATUS_BELONGSTOSET_GF   =  1,     /** all GUB variables are in noncovervars F (and noncovervars R) */
   GUBCONSSTATUS_BELONGSTOSET_GC2  =  2,     /** all GUB variables are in covervars C2 */
   GUBCONSSTATUS_BELONGSTOSET_GNC1 =  3,     /** some GUB variables are in covervars C1, others in noncovervars R or F */
   GUBCONSSTATUS_BELONGSTOSET_GOC1 =  4      /** all GUB variables are in covervars C1 */
};
typedef enum GUBConsstatus GUBCONSSTATUS;

/** data structure of GUB constraints */
struct SCIP_GUBCons
{
   int*                  gubvars;            /**< indices of GUB variables in knapsack constraint */
   GUBVARSTATUS*         gubvarsstatus;      /**< status of GUB variables */
   int                   ngubvars;           /**< number of GUB variables */
   int                   gubvarssize;        /**< size of gubvars array */
};
typedef struct SCIP_GUBCons SCIP_GUBCONS;

/** data structure of a set of GUB constraints */
struct SCIP_GUBSet
{
   SCIP_GUBCONS**        gubconss;           /**< GUB constraints in GUB set */
   GUBCONSSTATUS*        gubconsstatus;      /**< status of GUB constraints */
   int                   ngubconss;          /**< number of GUB constraints */
   int                   nvars;              /**< number of variables in knapsack constraint */
   int*                  gubconssidx;        /**< index of GUB constraint (in gubconss array) of each knapsack variable */
   int*                  gubvarsidx;         /**< index in GUB constraint (in gubvars array) of each knapsack variable  */
};
typedef struct SCIP_GUBSet SCIP_GUBSET;

/*
 * Local methods
 */

/** comparison method for two sorting key pairs */
static
SCIP_DECL_SORTPTRCOMP(compSortkeypairs)
{
   SORTKEYPAIR* sortkeypair1 = (SORTKEYPAIR*)elem1;
   SORTKEYPAIR* sortkeypair2 = (SORTKEYPAIR*)elem2;

   if( sortkeypair1->key1 < sortkeypair2->key1 )
      return -1;
   else if( sortkeypair1->key1 > sortkeypair2->key1 )
      return +1;
   else if( sortkeypair1->key2 < sortkeypair2->key2 )
      return -1;
   else if( sortkeypair1->key2 > sortkeypair2->key2 )
      return +1;
   else 
      return 0;
}

/** creates event data */
static
SCIP_RETCODE eventdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTDATA**      eventdata,          /**< pointer to store event data */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Longint          weight              /**< weight of variable */
   )
{
   assert(eventdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, eventdata) );
   (*eventdata)->cons = cons;
   (*eventdata)->weight = weight;

   return SCIP_OKAY;
}  

/** frees event data */
static
SCIP_RETCODE eventdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTDATA**      eventdata           /**< pointer to event data */
   )
{
   assert(eventdata != NULL);

   SCIPfreeBlockMemory(scip, eventdata);

   return SCIP_OKAY;
}

/** sorts items in knapsack with nonincreasing weights */
static
void sortItems(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->weights != NULL);
   assert(consdata->nvars == 0 || consdata->eventdata != NULL);
   assert(consdata->nvars == 0 || (consdata->cliquepartition != NULL && consdata->negcliquepartition != NULL));

   if( !consdata->sorted )
   {
      int pos;
      int lastcliquenum;
      int v;

      /* sort of five joint arrays of Long/pointer/pointer/ints/ints,
       * sorted by first array in non-increasing order via sort template */
      SCIPsortDownLongPtrPtrIntInt(
         consdata->weights,
         (void**)consdata->vars,
         (void**)consdata->eventdata,
         consdata->cliquepartition,
         consdata->negcliquepartition,
         consdata->nvars);

      v = consdata->nvars - 1;
      /* sort all items with same weight according to their variable index, used for hash value for fast pairwise comparison of all constraints */
      while( v >= 0 )
      {
         int w = v - 1;

         while( w >= 0 && consdata->weights[v] == consdata->weights[w] )
            --w;

         if( v - w > 1 )
         {
            /* sort all corresponding parts of arrays for which the weights are equal by using the variable index */
            SCIPsortPtrPtrIntInt(
               (void**)(&(consdata->vars[w+1])),
               (void**)(&(consdata->eventdata[w+1])),
               &(consdata->cliquepartition[w+1]),
               &(consdata->negcliquepartition[w+1]),
               SCIPvarComp,
               v - w);
         }
         v = w;
      }

      /* we need to make sure that our clique numbers of our normal clique will be in increasing order without gaps */
      if( consdata->cliquepartitioned )
      {
         lastcliquenum = 0;

         for( pos = 0; pos < consdata->nvars; ++pos )
         {
            /* if the clique number in the normal clique at position pos is greater than the last found clique number the
             * partition is invalid */
            if( consdata->cliquepartition[pos] > lastcliquenum )
            {
               consdata->cliquepartitioned = FALSE;
               break;
            }
            else if( consdata->cliquepartition[pos] == lastcliquenum )
               ++lastcliquenum;
         }
      }
      /* we need to make sure that our clique numbers of our negated clique will be in increasing order without gaps */
      if( consdata->negcliquepartitioned )
      {
         lastcliquenum = 0;

         for( pos = 0; pos < consdata->nvars; ++pos )
         {
            /* if the clique number in the negated clique at position pos is greater than the last found clique number the
             * partition is invalid */
            if( consdata->negcliquepartition[pos] > lastcliquenum )
            {
               consdata->negcliquepartitioned = FALSE;
               break;
            }
            else if( consdata->negcliquepartition[pos] == lastcliquenum )
               ++lastcliquenum;
         }
      }

      consdata->sorted = TRUE;
   }
#ifndef NDEBUG
   {
      /* check if the weight array is sorted in a non-increasing way */
      int i;
      for( i = 0; i < consdata->nvars-1; ++i )
         assert(consdata->weights[i] >= consdata->weights[i+1]);
   }
#endif
}

/** calculates a partition of the variables into cliques */
static
SCIP_RETCODE calcCliquepartition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< knapsack constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool             normalclique,       /**< Should normal cliquepartition be created? */
   SCIP_Bool             negatedclique       /**< Should negated cliquepartition be created? */
   )
{
   SCIP_Bool ispartitionoutdated;
   SCIP_Bool isnegpartitionoutdated;
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->cliquepartition != NULL && consdata->negcliquepartition != NULL));

   /* rerun eventually if number of global cliques increased considerably since last partition */
   ispartitionoutdated = (conshdlrdata->updatecliquepartitions && consdata->ncliques > 1
         && SCIPgetNCliques(scip) >= (int)(conshdlrdata->clqpartupdatefac * consdata->ncliqueslastpart));

   if( normalclique && ( !consdata->cliquepartitioned || ispartitionoutdated ) )
   {
      SCIP_CALL( SCIPcalcCliquePartition(scip, consdata->vars, consdata->nvars, consdata->cliquepartition, &consdata->ncliques) );
      consdata->cliquepartitioned = TRUE;
      consdata->ncliqueslastpart = SCIPgetNCliques(scip);
   }

   /* rerun eventually if number of global cliques increased considerably since last negated partition */
   isnegpartitionoutdated = (conshdlrdata->updatecliquepartitions && consdata->nnegcliques > 1
         && SCIPgetNCliques(scip) >= (int)(conshdlrdata->clqpartupdatefac * consdata->ncliqueslastnegpart));

   if( negatedclique && (!consdata->negcliquepartitioned || isnegpartitionoutdated) )
   {
      SCIP_CALL( SCIPcalcNegatedCliquePartition(scip, consdata->vars, consdata->nvars, consdata->negcliquepartition, &consdata->nnegcliques) );
      consdata->negcliquepartitioned = TRUE;
      consdata->ncliqueslastnegpart = SCIPgetNCliques(scip);
   }
   assert(!consdata->cliquepartitioned || consdata->ncliques <= consdata->nvars);
   assert(!consdata->negcliquepartitioned || consdata->nnegcliques <= consdata->nvars);

   return SCIP_OKAY;
}

/** installs rounding locks for the given variable in the given knapsack constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding up may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given knapsack constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding up may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** catches bound change events for variables in knapsack */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(cons != NULL);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->weights != NULL);
   assert(consdata->nvars == 0 || consdata->eventdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( eventdataCreate(scip, &consdata->eventdata[i], cons, consdata->weights[i]) );
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[i], EVENTTYPE_KNAPSACK,
            eventhdlr, consdata->eventdata[i], &consdata->eventdata[i]->filterpos) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for variables in knapsack */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->weights != NULL);
   assert(consdata->nvars == 0 || consdata->eventdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[i], EVENTTYPE_KNAPSACK,
            eventhdlr, consdata->eventdata[i], consdata->eventdata[i]->filterpos) );
      SCIP_CALL( eventdataFree(scip, &consdata->eventdata[i]) );
   }

   return SCIP_OKAY;
}

/** ensures, that vars and vals arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< knapsack constraint data */
   int                   num,                /**< minimum number of entries to store */
   SCIP_Bool             transformed         /**< is constraint from transformed problem? */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);

   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->weights, consdata->varssize, newsize) );
      if( transformed )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdata, consdata->varssize, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->cliquepartition, consdata->varssize, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->negcliquepartition, consdata->varssize, newsize) );
      }
      else
      {
         assert(consdata->eventdata == NULL);
         assert(consdata->cliquepartition == NULL);
         assert(consdata->negcliquepartition == NULL);
      }
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}

/** updates all weight sums for fixed and unfixed variables */
static
void updateWeightSums(
   SCIP_CONSDATA*        consdata,           /**< knapsack constraint data */
   SCIP_VAR*             var,                /**< variable for this weight */
   SCIP_Longint          weightdelta         /**< difference between the old and the new weight of the variable */
   )
{
   assert(consdata != NULL);
   assert(var != NULL);

   consdata->weightsum += weightdelta;

   if( SCIPvarGetLbLocal(var) > 0.5 )
      consdata->onesweightsum += weightdelta;

   assert(consdata->weightsum >= 0);
   assert(consdata->onesweightsum >= 0);
}

/** creates knapsack constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   int                   nvars,              /**< number of variables in knapsack */
   SCIP_VAR**            vars,               /**< variables of knapsack */
   SCIP_Longint*         weights,            /**< weights of knapsack items */
   SCIP_Longint          capacity            /**< capacity of knapsack */
   )
{
   int v;
   SCIP_Longint constant;

   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   constant = 0L;
   (*consdata)->vars = NULL;
   (*consdata)->weights = NULL;
   (*consdata)->nvars = 0;
   if( nvars > 0 )
   {
      SCIP_VAR** varsbuffer;
      SCIP_Longint* weightsbuffer;
      int k;

      SCIP_CALL( SCIPallocBufferArray(scip, &varsbuffer, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &weightsbuffer, nvars) );

      k = 0;
      for( v = 0; v < nvars; ++v )
      {
         assert(vars[v] != NULL);
         assert(SCIPvarIsBinary(vars[v]));

         /* all weight have to be non negative */
         assert( weights[v] >= 0 );

         if( weights[v] > 0 )
         {
            /* treat fixed variables as constants if problem compression is enabled */
            if( SCIPisConsCompressionEnabled(scip) && SCIPvarGetLbGlobal(vars[v]) > SCIPvarGetUbGlobal(vars[v]) - 0.5 )
            {
               /* only if the variable is fixed to 1, we add its weight to the constant */
               if( SCIPvarGetUbGlobal(vars[v]) > 0.5 )
                  constant += weights[v];
            }
            else
            {
               varsbuffer[k] = vars[v];
               weightsbuffer[k] = weights[v];
               ++k;
            }
         }
      }
      assert(k >= 0);

      (*consdata)->nvars = k;

      /* copy the active variables and weights into the constraint data structure */
      if( k > 0 )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, varsbuffer, k) );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->weights, weightsbuffer, k) );
      }

      /* free buffer storage */
      SCIPfreeBufferArray(scip, &weightsbuffer);
      SCIPfreeBufferArray(scip, &varsbuffer);
   }

   /* capacity has to be greater or equal to zero */
   assert(capacity >= 0);
   assert(constant >= 0);

   (*consdata)->varssize = (*consdata)->nvars;
   (*consdata)->capacity = capacity - constant;
   (*consdata)->eventdata = NULL;
   (*consdata)->cliquepartition = NULL;
   (*consdata)->negcliquepartition = NULL;
   (*consdata)->row = NULL;
   (*consdata)->weightsum = 0;
   (*consdata)->onesweightsum = 0;
   (*consdata)->ncliques = 0;
   (*consdata)->nnegcliques = 0;
   (*consdata)->presolvedtiming = 0;
   (*consdata)->sorted = FALSE;
   (*consdata)->cliquepartitioned = FALSE;
   (*consdata)->negcliquepartitioned = FALSE;
   (*consdata)->ncliqueslastpart = -1;
   (*consdata)->ncliqueslastnegpart = -1;
   (*consdata)->merged = FALSE;
   (*consdata)->cliquesadded = FALSE;
   (*consdata)->varsdeleted = FALSE;
   (*consdata)->existmultaggr = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      for( v = 0; v < (*consdata)->nvars; v++ )
      {
         SCIP_VAR* var = SCIPvarGetProbvar((*consdata)->vars[v]);
         assert(var != NULL);
         (*consdata)->existmultaggr = (*consdata)->existmultaggr || (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
      }

      /* allocate memory for additional data structures */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdata, (*consdata)->nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->cliquepartition, (*consdata)->nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->negcliquepartition, (*consdata)->nvars) );
   }

   /* calculate sum of weights and capture variables */
   for( v = 0; v < (*consdata)->nvars; ++v )
   {
      /* calculate sum of weights */
      updateWeightSums(*consdata, (*consdata)->vars[v], (*consdata)->weights[v]);

      /* capture variables */
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
   }
   return SCIP_OKAY;
}

/** frees knapsack constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to the constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }
   if( (*consdata)->eventdata != NULL )
   {
      SCIP_CALL( dropEvents(scip, *consdata, eventhdlr) );
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdata, (*consdata)->varssize);
   }
   if( (*consdata)->negcliquepartition != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->negcliquepartition, (*consdata)->varssize);
   }
   if( (*consdata)->cliquepartition != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->cliquepartition, (*consdata)->varssize);
   }
   if( (*consdata)->vars != NULL )
   {
      int v;

      /* release variables */
      for( v = 0; v < (*consdata)->nvars; v++ )
      {
         assert((*consdata)->vars[v] != NULL);
         SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
      }

      assert( (*consdata)->weights != NULL );
      assert( (*consdata)->varssize > 0 );
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->varssize);
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** changes a single weight in knapsack constraint data */
static
void consdataChgWeight(
   SCIP_CONSDATA*        consdata,           /**< knapsack constraint data */
   int                   item,               /**< item number */
   SCIP_Longint          newweight           /**< new weight of item */
   )
{
   SCIP_Longint oldweight;
   SCIP_Longint weightdiff;

   assert(consdata != NULL);
   assert(0 <= item && item < consdata->nvars);

   oldweight = consdata->weights[item];
   weightdiff = newweight - oldweight;
   consdata->weights[item] = newweight;


   /* update weight sums for all and fixed variables */
   updateWeightSums(consdata, consdata->vars[item], weightdiff);

   if( consdata->eventdata != NULL )
   {
      assert(consdata->eventdata[item] != NULL);
      assert(consdata->eventdata[item]->weight == oldweight);
      consdata->eventdata[item]->weight = newweight;
   }

   consdata->presolvedtiming = 0;
   consdata->sorted = FALSE;

   /* recalculate cliques extraction after a weight was increased */
   if( oldweight < newweight )
   {
      consdata->cliquesadded = FALSE;
   }
}

/** creates LP row corresponding to knapsack constraint */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< knapsack constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row, SCIPconsGetHdlr(cons), SCIPconsGetName(cons),
         -SCIPinfinity(scip), (SCIP_Real)consdata->capacity,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, consdata->row) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->vars[i], (SCIP_Real)consdata->weights[i]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, consdata->row) );

   return SCIP_OKAY;
}

/** adds linear relaxation of knapsack constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff              /**< whether a cutoff has been detected */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cutoff != NULL );
   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMsg(scip, "adding relaxation of knapsack constraint <%s> (capacity %" SCIP_LONGINT_FORMAT "): ",
         SCIPconsGetName(cons), consdata->capacity);
      SCIPdebug( SCIP_CALL(SCIPprintRow(scip, consdata->row, NULL)) );
      SCIP_CALL( SCIPaddRow(scip, consdata->row, FALSE, cutoff) );
   }

   return SCIP_OKAY;
}

/** checks knapsack constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "checking knapsack constraint <%s> for feasibility of solution %p (lprows=%u)\n",
      SCIPconsGetName(cons), (void*)sol, checklprows);

   *violated = FALSE;

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      SCIP_Real sum;
      SCIP_Longint integralsum;
      SCIP_Bool ishuge;
      SCIP_Real absviol;
      SCIP_Real relviol;
      int v;

      /* increase age of constraint; age is reset to zero, if a violation was found only in case we are in
       * enforcement
       */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }

      sum = 0.0;
      integralsum = 0;
      /* we perform a more exact comparison if the capacity does not exceed the huge value */
      if( SCIPisHugeValue(scip, (SCIP_Real) consdata->capacity) )
      {
         ishuge = TRUE;

         /* sum over all weight times the corresponding solution value */
         for( v = consdata->nvars - 1; v >= 0; --v )
         {
            assert(SCIPvarIsBinary(consdata->vars[v]));
            sum += consdata->weights[v] * SCIPgetSolVal(scip, sol, consdata->vars[v]);
         }
      }
      else
      {
         ishuge = FALSE;

         /* sum over all weight for which the variable has a solution value of 1 in feastol */
         for( v = consdata->nvars - 1; v >= 0; --v )
         {
            assert(SCIPvarIsBinary(consdata->vars[v]));

            if( SCIPgetSolVal(scip, sol, consdata->vars[v]) > 0.5 )
               integralsum += consdata->weights[v];
         }
      }

      /* calculate constraint violation and update it in solution */
      absviol = ishuge ? sum : (SCIP_Real)integralsum;
      absviol -= consdata->capacity;
      relviol = SCIPrelDiff(absviol + consdata->capacity, (SCIP_Real)consdata->capacity);
      if( sol != NULL )
         SCIPupdateSolLPConsViolation(scip, sol, absviol, relviol);

      if( SCIPisFeasPositive(scip, absviol) )
      {
         *violated = TRUE;

         /* only reset constraint age if we are in enforcement */
         if( sol == NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }

         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );

            SCIPinfoMessage(scip, NULL, ";\n");
            SCIPinfoMessage(scip, NULL, "violation: the capacity is violated by %.15g\n", absviol);
         }
      }

   }

   return SCIP_OKAY;
}

/* IDX computes the integer index for the optimal solution array */
#define IDX(j,d) ((j)*(intcap)+(d))

/** solves knapsack problem in maximization form exactly using dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 *
 * @note in case you provide the solitems or nonsolitems array you also have to provide the counter part as well
 */
SCIP_RETCODE SCIPsolveKnapsackExactly(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval,             /**< pointer to store optimal solution value, or NULL */
   SCIP_Bool*            success             /**< pointer to store if an error occured during solving
                                              *   (normally a memory problem) */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real* tempsort;
   SCIP_Real* optvalues;
   int intcap;
   int d;
   int j;
   SCIP_Longint weightsum;
   int* myitems;
   SCIP_Longint* myweights;
   int* allcurrminweight;
   SCIP_Real* myprofits;
   int nmyitems;
   SCIP_Longint gcd;
   SCIP_Longint minweight;
   SCIP_Longint maxweight;
   int currminweight;
   SCIP_Longint greedycap;
   SCIP_Longint greedysolweight;
   SCIP_Real greedysolvalue;
   SCIP_Bool eqweights;
   SCIP_Bool isoptimal;
   const size_t maxsize_t = (size_t)(-1);

   assert(weights != NULL);
   assert(profits != NULL);
   assert(capacity >= 0);
   assert(items != NULL);
   assert(nitems >= 0);
   assert(success != NULL);

   *success = TRUE;

#ifndef NDEBUG
   for( j = nitems - 1; j >= 0; --j )
      assert(weights[j] >= 0);
#endif

   SCIPdebugMsg(scip, "Solving knapsack exactly.\n");

   /* initializing solution value */
   if( solval != NULL )
      *solval = 0.0;

   /* produces optimal solution by following the table */
   if( solitems != NULL)
   {
      assert(items != NULL);
      assert(nsolitems != NULL);
      assert(nonsolitems != NULL);
      assert(nnonsolitems != NULL);

      *nnonsolitems = 0;
      *nsolitems = 0;
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &myweights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &myprofits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &myitems, nitems) );
   nmyitems = 0;
   weightsum = 0;
   minweight = SCIP_LONGINT_MAX;
   maxweight = 0;

   /* remove unnecessary items */
   for( j = 0; j < nitems; ++j )
   {
      assert(0 <= weights[j] && weights[j] < SCIP_LONGINT_MAX);
      /* items does not fit */
      if( weights[j] > capacity )
      {
         if( solitems != NULL)
         {
            nonsolitems[*nnonsolitems] = items[j];
            ++(*nnonsolitems);
         }
      }
      /* items we does not want */
      else if( profits[j] <= 0.0 )
      {
         if( solitems != NULL)
         {
            nonsolitems[*nnonsolitems] = items[j];
            ++(*nnonsolitems);
         }
      }
      /* items which always fit */
      else if( weights[j] == 0 )
      {
         if( solitems != NULL)
         {
            solitems[*nsolitems] = items[j];
            ++(*nsolitems);
         }
         if( solval != NULL )
            *solval += profits[j];
      }
      /* all important items */
      else
      {
         myweights[nmyitems] = weights[j];
         myprofits[nmyitems] = profits[j];
         myitems[nmyitems] = items[j];

         /* remember smallest item */
         if( myweights[nmyitems] < minweight )
            minweight = myweights[nmyitems];

         /* remember bigest item */
         if( myweights[nmyitems] > maxweight )
            maxweight = myweights[nmyitems];

         weightsum += myweights[nmyitems];
         ++nmyitems;
      }
   }

   /* no item is left then goto end */
   if( nmyitems == 0 )
   {
      SCIPdebugMsg(scip, "After preprocessing no items are left.\n");

      goto TERMINATE;
   }
   /* if all items fit, we also do not need to do the expensive stuff later on */
   else if( weightsum > 0 && weightsum <= capacity )
   {
      SCIPdebugMsg(scip, "After preprocessing all items fit into knapsack.\n");

      for( j = nmyitems - 1; j >= 0; --j )
      {
         if( solitems != NULL )
         {
            solitems[*nsolitems] = myitems[j];
            ++(*nsolitems);
         }
         if( solval != NULL )
            *solval += myprofits[j];
      }

      goto TERMINATE;
   }

   assert(minweight > 0);
   assert(maxweight > 0);

   if( maxweight > 1 )
   {
      /* determine greatest common divisor */
      gcd = myweights[nmyitems - 1];
      for( j = nmyitems - 2; j >= 0 && gcd >= 2; --j )
         gcd = SCIPcalcGreComDiv(gcd, myweights[j]);

      SCIPdebugMsg(scip, "Gcd is %" SCIP_LONGINT_FORMAT ".\n", gcd);

      /* divide by greatest common divisor */
      if( gcd > 1 )
      {
         eqweights = TRUE;
         for( j = nmyitems - 1; j >= 0; --j )
         {
            myweights[j] /= gcd;
            eqweights = eqweights && (myweights[j] == 1);
         }
         capacity /= gcd;
         minweight /= gcd;
      }
      else
         eqweights = FALSE;
   }
   else
   {
      assert(maxweight == 1);
      eqweights = TRUE;
   }

   assert(minweight <= capacity);

   /* only one item fits, than take the best */
   if( minweight > capacity / 2 )
   {
      int p;

      SCIPdebugMsg(scip, "Only one item fits into knapsack, so take the best.\n");

      p = nmyitems - 1;

      /* find best item */
      for( j = nmyitems - 2; j >= 0; --j )
         if( myprofits[j] > myprofits[p] )
            p = j;

      /* update solution information */
      if( solitems != NULL)
      {
         solitems[*nsolitems] = myitems[p];
         ++(*nsolitems);
         for( j = nmyitems - 1; j >= 0; --j )
            if( j != p )
            {
               nonsolitems[*nnonsolitems] = myitems[j];
               ++(*nnonsolitems);
            }
      }
      /* update solution value */
      if( solval != NULL )
         *solval += myprofits[p];

      goto TERMINATE;
   }

   /* all items have the same weight, than take the best */
   if( eqweights )
   {
      SCIP_Real addval;

      SCIPdebugMsg(scip, "All weights are equal, so take the best.\n");

      SCIPsortDownRealIntLong(myprofits, myitems, myweights, nmyitems);

      addval = 0.0;
      /* update solution information */
      if( solitems != NULL || solval != NULL )
      {
         SCIP_Longint i;

         /* if all items would fit we had handled this case before */
         assert((SCIP_Longint) nmyitems > capacity);

         /* take the first best items into the solution */
         for( i = capacity - 1; i >= 0; --i )
         {
            if( solitems != NULL)
            {
               assert(nonsolitems != NULL);
               solitems[*nsolitems] = myitems[i];
               ++(*nsolitems);
            }
            addval += myprofits[i];
         }

         if( solitems != NULL)
         {
            assert(nonsolitems != NULL);

            /* the rest are not in the solution */
            for( i = nmyitems - 1; i >= capacity; --i )
            {
               nonsolitems[*nnonsolitems] = myitems[i];
               ++(*nnonsolitems);
            }
         }
      }
      /* update solution value */
      if( solval != NULL )
      {
         assert(addval > 0.0);
         *solval += addval;
      }

      goto TERMINATE;
   }

   /* in the following table we do not need the first minweight columns */
   capacity -= (minweight - 1);

   /* we can only handle integers */
   if( capacity >= INT_MAX )
   {
      SCIPdebugMsg(scip, "Capacity is to big, so we cannot handle it here.\n");

      *success = FALSE;
      goto TERMINATE;
   }
   assert(capacity < INT_MAX);

   intcap = (int)capacity;
   assert(intcap >= 0);
   assert(nmyitems > 0);
   assert(sizeof(size_t) >= sizeof(int)); /* no following conversion should be messed up */

   /* this condition checks if we will try to allocate a correct number of bytes and do not have an overflow, while
    * computing the size for the allocation
    */
   if( intcap < 0 || (intcap > 0 && (((size_t)nmyitems) > (maxsize_t / (size_t)intcap / sizeof(*optvalues)) || ((size_t)nmyitems) * ((size_t)intcap) * sizeof(*optvalues) > ((size_t)INT_MAX) )) ) /*lint !e571*/
   {
      SCIPdebugMsg(scip, "Too much memory (%lu) would be consumed.\n", (unsigned long) (((size_t)nmyitems) * ((size_t)intcap) * sizeof(*optvalues))); /*lint !e571*/

      *success = FALSE;
      goto TERMINATE;
   }

   /* allocate temporary memory and check for memory exceeding */
   retcode = SCIPallocBufferArray(scip, &optvalues, nmyitems * intcap);
   if( retcode == SCIP_NOMEMORY )
   {
      SCIPdebugMsg(scip, "Did not get enough memory.\n");

      *success = FALSE;
      goto TERMINATE;
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* sort myitems (plus corresponding arrays myweights and myprofits) such that
    * p_1/w_1 >= p_2/w_2 >= ... >= p_n/w_n, this is only use for greedy solution
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nmyitems) );
   for( j = nmyitems - 1; j >= 0; --j )
      tempsort[j] = myprofits[j]/((SCIP_Real) myweights[j]);

   SCIPsortDownRealLongRealInt(tempsort, myweights, myprofits, myitems, nmyitems);

   /* initialize values for greedy solution information */
   greedysolweight = 0;
   greedysolvalue = 0.0;
   isoptimal = TRUE;
   greedycap = capacity + (minweight - 1);

   SCIPdebugMsg(scip, "Determine greedy solution.\n");

   /* determine greedy solution */
   for( j = 0; j < nmyitems; ++j )
   {
      assert(myweights[j] <= greedycap);

      /* take all fitting items */
      if( myweights[j] + greedysolweight <= greedycap )
      {
         /* update greedy solution weight and value */
         greedysolweight += myweights[j];
         greedysolvalue += myprofits[j];
         continue;
      }
      else if( greedysolweight < greedycap )
         isoptimal = FALSE;
      break;
   }
   assert(greedysolweight > 0);
   assert(greedysolvalue > 0.0);

   /* greedy solution is optimal */
   if( isoptimal )
   {
      assert(greedysolweight == greedycap);

      SCIPdebugMsg(scip, "Greedy solution is optimal.\n");

      greedysolweight = 0;

      /* update solution information */
      if( solitems != NULL)
      {
         /* take the first best items into the solution */
         for( j = 0; j < nmyitems; ++j )
         {
            /* take all fitting items */
            if( myweights[j] + greedysolweight <= greedycap )
            {
               solitems[*nsolitems] = myitems[j];
               ++(*nsolitems);
               greedysolweight += myweights[j];
            }
            else
            {
               nonsolitems[*nnonsolitems] = myitems[j];
               ++(*nnonsolitems);
            }
         }
      }
      /* update solution value */
      if( solval != NULL )
      {
         assert(greedysolvalue > 0.0);
         *solval += greedysolvalue;
      }

      SCIPfreeBufferArray(scip, &tempsort);
      SCIPfreeBufferArray(scip, &optvalues);

      goto TERMINATE;
   }

   SCIPdebugMsg(scip, "Start real exact algorithm.\n");

   /* we memorize at each step the current minimal weight to later on know which value in our optvalues matrix is valid;
    * all values entries of the j-th row of optvalues is valid if the index is >= allcurrminweight[j], otherwise it is
    * invalid, a second possibility would be to clear the whole optvalues, which should be more expensive than storing
    * 'nmyitem' values
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &allcurrminweight, nmyitems) );
   assert(myweights[0] - minweight < INT_MAX);
   currminweight = (int) (myweights[0] - minweight);
   allcurrminweight[0] = currminweight;

   /* fills first row of dynamic programming table with optimal values */
   for( d = currminweight; d < intcap; ++d )
      optvalues[d] = myprofits[0];
   /* fills dynamic programming table with optimal values */
   for( j = 1; j < nmyitems; ++j )
   {
      int intweight;

      /* compute important part of weight, which will be represented in the table */
      intweight = (int)(myweights[j] - minweight);
      assert(0 <= intweight && intweight < intcap);

      /* copy all nonzeros from row above */
      for( d = currminweight; d < intweight && d < intcap; ++d )
         optvalues[IDX(j,d)] = optvalues[IDX(j-1,d)];

      /* update corresponding row */
      for( d = intweight; d < intcap; ++d )
      {
         /* if index d is smaller the the current minweight then optvalues[IDX(j-1,d)] is not initialized, i.e. should
          * be 0
          */
         if( d < currminweight )
         {
            optvalues[IDX(j,d)] = myprofits[j];
         }
         else
         {
            SCIP_Real sumprofit;

            if( d - myweights[j] < currminweight )
               sumprofit = myprofits[j];
            else
               sumprofit = optvalues[IDX(j-1,(int)(d-myweights[j]))] + myprofits[j];

            optvalues[IDX(j,d)] = MAX(sumprofit, optvalues[IDX(j-1,d)]);
         }
      }
      /* update currminweight */
      if( intweight < currminweight )
         currminweight = intweight;

      allcurrminweight[j] = currminweight;
   }

   /* update optimal solution by following the table */
   if( solitems != NULL)
   {
      d = intcap - 1;

      SCIPdebugMsg(scip, "Fill the solution vector after solving exactly.\n");

      /* insert all items in (non-) solution vector */
      for( j = nmyitems - 1; j > 0; --j )
      {
         /* if we cannot find any item anymore which is in our solution stop, if the following condition holds this
          * means all remaining items does not fit anymore
          */
         if( d < allcurrminweight[j] )
         {
            /* we cannot have exceeded our capacity */
            assert((SCIP_Longint) d >= -minweight);
            break;
         }
         /* collect solution items, first condition means that no next item can fit anymore, but this does */
         if( d < allcurrminweight[j-1] || optvalues[IDX(j,d)] > optvalues[IDX(j-1,d)] )
         {
            solitems[*nsolitems] = myitems[j];
            ++(*nsolitems);

            /* check that we do not have an underflow */
            assert(myweights[j] <= (INT_MAX + (SCIP_Longint) d));
            d = (int)(d - myweights[j]);
         }
         /* collect non-solution items */
         else
         {
            nonsolitems[*nnonsolitems] = myitems[j];
            ++(*nnonsolitems);
         }
      }

      /* insert remaining items */
      if( d >= allcurrminweight[j] )
      {
         assert(j == 0);
         solitems[*nsolitems] = myitems[j];
         ++(*nsolitems);
      }
      else
      {
         assert(j >= 0);
         assert(d < allcurrminweight[j]);

         for( ; j >= 0; --j )
         {
            nonsolitems[*nnonsolitems] = myitems[j];
            ++(*nnonsolitems);
         }
      }

      assert(*nsolitems + *nnonsolitems == nitems);
   }

   /* update solution value */
   if( solval != NULL )
      *solval += optvalues[IDX(nmyitems-1,intcap-1)];

   SCIPfreeBufferArray(scip, &allcurrminweight);

   /* free all temporary memory */
   SCIPfreeBufferArray(scip, &tempsort);
   SCIPfreeBufferArray(scip, &optvalues);

 TERMINATE:
   SCIPfreeBufferArray(scip, &myitems);
   SCIPfreeBufferArray(scip, &myprofits);
   SCIPfreeBufferArray(scip, &myweights);

   return SCIP_OKAY;
}

/** solves knapsack problem in maximization form approximately by solving the LP-relaxation of the problem using Dantzig's
 *  method and rounding down the solution; if needed, one can provide arrays to store all selected items and all not
 *  selected items
 */
SCIP_RETCODE SCIPsolveKnapsackApproximately(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   )
{
   SCIP_Real* tempsort;
   SCIP_Longint solitemsweight;
   SCIP_Real* realweights;
   int j;
   int criticalindex;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(capacity >= 0);
   assert(items != NULL);
   assert(nitems >= 0);

   if( solitems != NULL )
   {
      *nsolitems = 0;
      *nnonsolitems = 0;
   }
   if( solval != NULL )
      *solval = 0.0;

   /* initialize data for median search */
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &realweights, nitems) );
   for( j = nitems - 1; j >= 0; --j )
   {
      tempsort[j] = profits[j]/((SCIP_Real) weights[j]);
      realweights[j] = (SCIP_Real)weights[j];

   }

   /* partially sort indices such that all elements that are larger than the break item appear first */
   SCIPselectWeightedDownRealLongRealInt(tempsort, weights, profits, items, realweights, (SCIP_Real)capacity, nitems, &criticalindex);

   /* selects items as long as they fit into the knapsack */
   solitemsweight = 0;
   for( j = 0; j < nitems && solitemsweight + weights[j] <= capacity; ++j )
   {
      if( solitems != NULL )
      {
         solitems[*nsolitems] = items[j];
         (*nsolitems)++;
      }
      if( solval != NULL )
         (*solval) += profits[j];
      solitemsweight += weights[j];
   }
   for( ; j < nitems && solitems != NULL; j++ )
   {
      nonsolitems[*nnonsolitems] = items[j];
      (*nnonsolitems)++;
   }

   SCIPfreeBufferArray(scip, &realweights);
   SCIPfreeBufferArray(scip, &tempsort);

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** prints all nontrivial GUB constraints and their LP solution values */
static
void GUBsetPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   SCIP_Real*            solvals             /**< solution values of variables in knapsack constraint; or NULL */
   )
{
   int nnontrivialgubconss;
   int c;

   nnontrivialgubconss = 0;

   SCIPdebugMsg(scip, "   Nontrivial GUBs of current GUB set:\n");

   /* print out all nontrivial GUB constraints, i.e., with more than one variable */
   for( c = 0; c < gubset->ngubconss; c++ )
   {
      SCIP_Real gubsolval;

      assert(gubset->gubconss[c]->ngubvars >= 0);

      /* nontrivial GUB */
      if( gubset->gubconss[c]->ngubvars > 1 )
      {
         int v;

         gubsolval = 0.0;
         SCIPdebugMsg(scip, "   GUB<%d>:\n", c);

         /* print GUB var */
         for( v = 0; v < gubset->gubconss[c]->ngubvars; v++ )
         {
            int currentvar;

            currentvar = gubset->gubconss[c]->gubvars[v];
            if( solvals != NULL )
            {
               gubsolval += solvals[currentvar];
               SCIPdebugMsg(scip, "      +<%s>(%4.2f)\n", SCIPvarGetName(vars[currentvar]), solvals[currentvar]);
            }
            else
            {
               SCIPdebugMsg(scip, "      +<%s>\n", SCIPvarGetName(vars[currentvar]));
            }
         }

	 /* check whether LP solution satisfies the GUB constraint */
         if( solvals != NULL )
         {
            SCIPdebugMsg(scip, "      =%4.2f <= 1 %s\n", gubsolval,
               SCIPisFeasGT(scip, gubsolval, 1.0) ? "--> violated" : "");
         }
         else
         {
            SCIPdebugMsg(scip, "      <= 1 %s\n", SCIPisFeasGT(scip, gubsolval, 1.0) ? "--> violated" : "");
         }
         nnontrivialgubconss++;
      }
   }

   SCIPdebugMsg(scip, "   --> %d/%d nontrivial GUBs\n", nnontrivialgubconss, gubset->ngubconss);
}
#endif

/** creates an empty GUB constraint */
static
SCIP_RETCODE GUBconsCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBCONS**        gubcons             /**< pointer to store GUB constraint data */
   )
{
   assert(scip != NULL);
   assert(gubcons != NULL);

   /* allocate memory for GUB constraint data structures */
   SCIP_CALL( SCIPallocBuffer(scip, gubcons) );
   (*gubcons)->gubvarssize = GUBCONSGROWVALUE;
   SCIP_CALL( SCIPallocBufferArray(scip, &(*gubcons)->gubvars, (*gubcons)->gubvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*gubcons)->gubvarsstatus, (*gubcons)->gubvarssize) );

   (*gubcons)->ngubvars = 0;

   return SCIP_OKAY;
}

/** frees GUB constraint */
static
SCIP_RETCODE GUBconsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBCONS**        gubcons             /**< pointer to GUB constraint data structure */
   )
{
   assert(scip != NULL);
   assert(gubcons != NULL);
   assert((*gubcons)->gubvars != NULL);
   assert((*gubcons)->gubvarsstatus != NULL);

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &(*gubcons)->gubvarsstatus);
   SCIPfreeBufferArray(scip, &(*gubcons)->gubvars);
   SCIPfreeBuffer(scip, gubcons);

   return SCIP_OKAY;
}

/** adds variable to given GUB constraint */
static
SCIP_RETCODE GUBconsAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBCONS*         gubcons,            /**< GUB constraint data */
   int                   var                 /**< index of given variable in knapsack constraint */
   )
{
   assert(scip != NULL);
   assert(gubcons != NULL);
   assert(gubcons->ngubvars >= 0 && gubcons->ngubvars < gubcons->gubvarssize);
   assert(gubcons->gubvars != NULL);
   assert(gubcons->gubvarsstatus != NULL);
   assert(var >= 0);

   /* add variable to GUB constraint */
   gubcons->gubvars[gubcons->ngubvars] = var;
   gubcons->gubvarsstatus[gubcons->ngubvars] = GUBVARSTATUS_UNINITIAL;
   gubcons->ngubvars++;

   /* increase space allocated to GUB constraint if the number of variables reaches the size */
   if( gubcons->ngubvars == gubcons->gubvarssize )
   {
      int newlen;

      newlen = gubcons->gubvarssize + GUBCONSGROWVALUE;
      SCIP_CALL( SCIPreallocBufferArray(scip, &gubcons->gubvars, newlen) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &gubcons->gubvarsstatus, newlen) );

      gubcons->gubvarssize = newlen;
   }

   return SCIP_OKAY;
}

/** deletes variable from its current GUB constraint */
static
SCIP_RETCODE GUBconsDelVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBCONS*         gubcons,            /**< GUB constraint data */
   int                   var,                /**< index of given variable in knapsack constraint */
   int                   gubvarsidx          /**< index of the variable in its current GUB constraint */
   )
{
   assert(scip != NULL);
   assert(gubcons != NULL);
   assert(var >= 0);
   assert(gubvarsidx >= 0 && gubvarsidx < gubcons->ngubvars);
   assert(gubcons->ngubvars >= gubvarsidx+1);
   assert(gubcons->gubvars[gubvarsidx] == var);

   /* delete variable from GUB by swapping it replacing in by the last variable in the GUB constraint */
   gubcons->gubvars[gubvarsidx] = gubcons->gubvars[gubcons->ngubvars-1];
   gubcons->gubvarsstatus[gubvarsidx] = gubcons->gubvarsstatus[gubcons->ngubvars-1];
   gubcons->ngubvars--;

   /* decrease space allocated for the GUB constraint, if the last GUBCONSGROWVALUE+1 array entries are now empty */
   if( gubcons->ngubvars < gubcons->gubvarssize - GUBCONSGROWVALUE && gubcons->ngubvars > 0 )
   {
      int newlen;

      newlen = gubcons->gubvarssize - GUBCONSGROWVALUE;

      SCIP_CALL( SCIPreallocBufferArray(scip, &gubcons->gubvars, newlen) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &gubcons->gubvarsstatus, newlen) );

      gubcons->gubvarssize = newlen;
   }

   return SCIP_OKAY;
}

/** moves variable from current GUB constraint to a different existing (nonempty) GUB constraint */
static
SCIP_RETCODE GUBsetMoveVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   var,                /**< index of given variable in knapsack constraint */
   int                   oldgubcons,         /**< index of old GUB constraint of given variable */
   int                   newgubcons          /**< index of new GUB constraint of given variable */
   )
{
   int oldgubvaridx;
   int replacevar;
   int j;

   assert(scip != NULL);
   assert(gubset != NULL);
   assert(var >= 0);
   assert(oldgubcons >= 0 && oldgubcons < gubset->ngubconss);
   assert(newgubcons >= 0 && newgubcons < gubset->ngubconss);
   assert(oldgubcons != newgubcons);
   assert(gubset->gubconssidx[var] == oldgubcons);
   assert(gubset->gubconss[oldgubcons]->ngubvars > 0);
   assert(gubset->gubconss[newgubcons]->ngubvars >= 0);

   SCIPdebugMsg(scip, "   moving variable<%s> from GUB<%d> to GUB<%d>\n", SCIPvarGetName(vars[var]), oldgubcons, newgubcons);

   oldgubvaridx = gubset->gubvarsidx[var];

   /* delete variable from old GUB constraint by replacing it by the last variable of the GUB constraint */
   SCIP_CALL( GUBconsDelVar(scip, gubset->gubconss[oldgubcons], var, oldgubvaridx) );

   /* in GUB set, update stored index of variable in old GUB constraint for the variable used for replacement;
    * replacement variable is given by old position of the deleted variable
    */
   replacevar = gubset->gubconss[oldgubcons]->gubvars[oldgubvaridx];
   assert(gubset->gubvarsidx[replacevar] == gubset->gubconss[oldgubcons]->ngubvars);
   gubset->gubvarsidx[replacevar] = oldgubvaridx;

   /* add variable to the end of new GUB constraint */
   SCIP_CALL( GUBconsAddVar(scip, gubset->gubconss[newgubcons], var) );
   assert(gubset->gubconss[newgubcons]->gubvars[gubset->gubconss[newgubcons]->ngubvars-1] == var);

   /* in GUB set, update stored index of GUB of moved variable and stored index of variable in this GUB constraint */
   gubset->gubconssidx[var] = newgubcons;
   gubset->gubvarsidx[var] = gubset->gubconss[newgubcons]->ngubvars-1;

   /* delete old GUB constraint if it became empty */
   if( gubset->gubconss[oldgubcons]->ngubvars == 0 )
   {
      SCIPdebugMsg(scip, "deleting empty GUB cons<%d> from current GUB set\n", oldgubcons);
#ifdef SCIP_DEBUG
      GUBsetPrint(scip, gubset, vars, NULL);
#endif

      /* free old GUB constraint */
      SCIP_CALL( GUBconsFree(scip, &gubset->gubconss[oldgubcons]) );

      /* if empty GUB was not the last one in GUB set data structure, replace it by last GUB constraint */
      if( oldgubcons != gubset->ngubconss-1 )
      {
         gubset->gubconss[oldgubcons] = gubset->gubconss[gubset->ngubconss-1];
         gubset->gubconsstatus[oldgubcons] = gubset->gubconsstatus[gubset->ngubconss-1];

         /* in GUB set, update stored index of GUB constraint for all variable of the GUB constraint used for replacement;
          * replacement GUB is given by old position of the deleted GUB
          */
         for( j = 0; j < gubset->gubconss[oldgubcons]->ngubvars; j++ )
         {
            assert(gubset->gubconssidx[gubset->gubconss[oldgubcons]->gubvars[j]] == gubset->ngubconss-1);
            gubset->gubconssidx[gubset->gubconss[oldgubcons]->gubvars[j]] = oldgubcons;
         }
      }

      /* update number of GUB constraints */
      gubset->ngubconss--;

      /* variable should be at given new position, unless new GUB constraint replaced empty old GUB constraint
       * (because it was at the end of the GUB constraint array)
       */
      assert(gubset->gubconssidx[var] == newgubcons
         || (newgubcons == gubset->ngubconss && gubset->gubconssidx[var] == oldgubcons));
   }
#ifndef NDEBUG
   else
      assert(gubset->gubconssidx[var] == newgubcons);
#endif

   return SCIP_OKAY;
}

/** swaps two variables in the same GUB constraint */
static
void GUBsetSwapVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure */
   int                   var1,               /**< first variable to be swapped */
   int                   var2                /**< second variable to be swapped */
   )
{
   int gubcons;
   int var1idx;
   GUBVARSTATUS var1status;
   int var2idx;
   GUBVARSTATUS var2status;

   assert(scip != NULL);
   assert(gubset != NULL);

   gubcons = gubset->gubconssidx[var1];
   assert(gubcons == gubset->gubconssidx[var2]);

   /* nothing to be done if both variables are the same */
   if( var1 == var2 )
      return;

   /* swap index and status of variables in GUB constraint */
   var1idx = gubset->gubvarsidx[var1];
   var1status = gubset->gubconss[gubcons]->gubvarsstatus[var1idx];
   var2idx = gubset->gubvarsidx[var2];
   var2status = gubset->gubconss[gubcons]->gubvarsstatus[var2idx];

   gubset->gubvarsidx[var1] = var2idx;
   gubset->gubconss[gubcons]->gubvars[var1idx] = var2;
   gubset->gubconss[gubcons]->gubvarsstatus[var1idx] = var2status;

   gubset->gubvarsidx[var2] = var1idx;
   gubset->gubconss[gubcons]->gubvars[var2idx] = var1;
   gubset->gubconss[gubcons]->gubvarsstatus[var2idx] = var1status;
}

/** initializes partition of knapsack variables into nonoverlapping trivial GUB constraints (GUB with one variable) */
static
SCIP_RETCODE GUBsetCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET**         gubset,             /**< pointer to store GUB set data structure */
   int                   nvars,              /**< number of variables in the knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity            /**< capacity of knapsack */
   )
{
   int i;

   assert(scip != NULL);
   assert(gubset != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(capacity >= 0);

   /* allocate memory for GUB set data structures */
   SCIP_CALL( SCIPallocBuffer(scip, gubset) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*gubset)->gubconss, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*gubset)->gubconsstatus, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*gubset)->gubconssidx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*gubset)->gubvarsidx, nvars) );
   (*gubset)->ngubconss = nvars;
   (*gubset)->nvars = nvars;

   /* initialize the set of GUB constraints */
   for( i = 0; i < nvars; i++ )
   {
      /* assign each variable to a new (trivial) GUB constraint */
      SCIP_CALL( GUBconsCreate(scip, &(*gubset)->gubconss[i]) );
      SCIP_CALL( GUBconsAddVar(scip, (*gubset)->gubconss[i], i) );

      /* set status of GUB constraint to initial */
      (*gubset)->gubconsstatus[i] = GUBCONSSTATUS_UNINITIAL;

      (*gubset)->gubconssidx[i] = i;
      (*gubset)->gubvarsidx[i] = 0;
      assert((*gubset)->gubconss[i]->ngubvars == 1);

      /* already updated status of variable in GUB constraint if it exceeds the capacity of the knapsack */
      if( weights[i] > capacity )
         (*gubset)->gubconss[(*gubset)->gubconssidx[i]]->gubvarsstatus[(*gubset)->gubvarsidx[i]] = GUBVARSTATUS_CAPACITYEXCEEDED;

   }

   return SCIP_OKAY;
}

/** frees GUB set data structure */
static
SCIP_RETCODE GUBsetFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET**         gubset              /**< pointer to GUB set data structure */
   )
{
   int i;

   assert(scip != NULL);
   assert(gubset != NULL);
   assert((*gubset)->gubconss != NULL);
   assert((*gubset)->gubconsstatus != NULL);
   assert((*gubset)->gubconssidx != NULL);
   assert((*gubset)->gubvarsidx != NULL);

   /* free all GUB constraints */
   for( i = (*gubset)->ngubconss-1; i >= 0; --i )
   {
      assert((*gubset)->gubconss[i] != NULL);
      SCIP_CALL( GUBconsFree(scip, &(*gubset)->gubconss[i]) );
   }

   /* free allocated memory */
   SCIPfreeBufferArray( scip, &(*gubset)->gubvarsidx );
   SCIPfreeBufferArray( scip, &(*gubset)->gubconssidx );
   SCIPfreeBufferArray( scip, &(*gubset)->gubconsstatus );
   SCIPfreeBufferArray( scip, &(*gubset)->gubconss );
   SCIPfreeBuffer(scip, gubset);

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks whether GUB set data structure is consistent */
static
SCIP_RETCODE GUBsetCheck(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure */
   SCIP_VAR**            vars                /**< variables in the knapsack constraint */
   )
{
   int i;
   int gubconsidx;
   int gubvaridx;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_Bool var1negated;
   SCIP_Bool var2negated;

   assert(scip != NULL);
   assert(gubset != NULL);

   SCIPdebugMsg(scip, "   GUB set consistency check:\n");

   /* checks for all knapsack vars consistency of stored index of associated gubcons and corresponding index in gubvars */
   for( i = 0; i < gubset->nvars; i++ )
   {
      gubconsidx = gubset->gubconssidx[i];
      gubvaridx = gubset->gubvarsidx[i];

      if( gubset->gubconss[gubconsidx]->gubvars[gubvaridx] != i )
      {
	 SCIPdebugMsg(scip, "   var<%d> should be in GUB<%d> at position<%d>, but stored is var<%d> instead\n", i,
            gubconsidx, gubvaridx, gubset->gubconss[gubconsidx]->gubvars[gubvaridx] );
      }
      assert(gubset->gubconss[gubconsidx]->gubvars[gubvaridx] == i);
   }

   /* checks for each GUB whether all pairs of its variables have a common clique */
   for( i = 0; i < gubset->ngubconss; i++ )
   {
      int j;

      for( j = 0; j < gubset->gubconss[i]->ngubvars; j++ )
      {
         int k;

         /* get corresponding active problem variable */
         var1 = vars[gubset->gubconss[i]->gubvars[j]];
         var1negated = FALSE;
         SCIP_CALL( SCIPvarGetProbvarBinary(&var1, &var1negated) );

         for( k = j+1; k < gubset->gubconss[i]->ngubvars; k++ )
         {
            /* get corresponding active problem variable */
            var2 = vars[gubset->gubconss[i]->gubvars[k]];
            var2negated = FALSE;
            SCIP_CALL( SCIPvarGetProbvarBinary(&var2, &var2negated) );

            if( !SCIPvarsHaveCommonClique(var1, !var1negated, var2, !var2negated, TRUE) )
            {
               SCIPdebugMsg(scip, "   GUB<%d>: var<%d,%s> and var<%d,%s> do not share a clique\n", i, j,
                  SCIPvarGetName(vars[gubset->gubconss[i]->gubvars[j]]), k,
                  SCIPvarGetName(vars[gubset->gubconss[i]->gubvars[k]]));
               SCIPdebugMsg(scip, "   GUB<%d>: var<%d,%s> and var<%d,%s> do not share a clique\n", i, j,
                  SCIPvarGetName(var1), k,
                  SCIPvarGetName(var2));
            }

            /* @todo: in case we used also negated cliques for the GUB partition, this assert has to be changed */
            assert(SCIPvarsHaveCommonClique(var1, !var1negated, var2, !var2negated, TRUE));
         }
      }
   }
   SCIPdebugMsg(scip, "   --> successful\n");

   return SCIP_OKAY;
}
#endif

/** calculates a partition of the given set of binary variables into cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same clique;
 *  the first variable is always assigned to clique 0, and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  note: in contrast to SCIPcalcCliquePartition(), variables with LP value 1 are put into trivial cliques (with one
 *  variable) and for the remaining variables, a partition with a small number of cliques is constructed
 */

static
SCIP_RETCODE GUBsetCalcCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques,           /**< pointer to store number of cliques actually contained in the partition */
   SCIP_Real*            solvals             /**< solution values of all given binary variables */
   )
{
   SCIP_VAR** tmpvars;
   SCIP_VAR** cliquevars;
   SCIP_Bool* cliquevalues;
   SCIP_Bool* tmpvalues;
   int* varseq;
   int* sortkeys;
   int ncliquevars;
   int maxncliquevarscomp;
   int nignorevars;
   int nvarsused;
   int i;

   assert(scip != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || cliquepartition != NULL);
   assert(ncliques != NULL);

   if( nvars == 0 )
   {
      *ncliques = 0;
      return SCIP_OKAY;
   }

   /* allocate temporary memory for storing the variables of the current clique */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquevalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvalues, nvars) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpvars, vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varseq, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeys, nvars) );

   /* initialize the cliquepartition array with -1 */
   /* initialize the tmpvalues array */
   for( i = nvars - 1; i >= 0; --i )
   {
      tmpvalues[i] = TRUE;
      cliquepartition[i] = -1;
   }

   /* get corresponding active problem variables */
   SCIP_CALL( SCIPvarsGetProbvarBinary(&tmpvars, &tmpvalues, nvars) );

   /* ignore variables with LP value 1 (will be assigned to trivial GUBs at the end) and sort remaining variables
    * by nondecreasing number of cliques the variables are in
    */
   nignorevars = 0;
   nvarsused = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( SCIPisFeasEQ(scip, solvals[i], 1.0) )
      {
         /* variables with LP value 1 are put to the end of varseq array and will not be sorted */
         varseq[nvars-1-nignorevars] = i;
         nignorevars++;
      }
      else
      {
         /* remaining variables are put to the front of varseq array and will be sorted by their number of cliques */
         varseq[nvarsused] = i;
         sortkeys[nvarsused] = SCIPvarGetNCliques(tmpvars[i], tmpvalues[i]);
         nvarsused++;
      }
   }
   assert(nvarsused + nignorevars == nvars);

   /* sort variables with LP value less than 1 by nondecreasing order of the number of cliques they are in */
   SCIPsortIntInt(sortkeys, varseq, nvarsused);

   maxncliquevarscomp = MIN(nvars*nvars, MAXNCLIQUEVARSCOMP);

   /* calculate the clique partition */
   *ncliques = 0;
   for( i = 0; i < nvars; ++i )
   {
      if( cliquepartition[varseq[i]] == -1 )
      {
         int j;

         /* variable starts a new clique */
         cliquepartition[varseq[i]] = *ncliques;
         cliquevars[0] = tmpvars[varseq[i]];
         cliquevalues[0] = tmpvalues[varseq[i]];
         ncliquevars = 1;

         /* if variable is not active (multi-aggregated or fixed), it cannot be in any clique and
          * if the variable has LP value 1 we do not want it to be in nontrivial cliques
          */
         if( SCIPvarIsActive(tmpvars[varseq[i]]) && i < nvarsused )
         {
            /* greedily fill up the clique */
            for( j = i + 1; j < nvarsused; ++j )
            {
               /* if variable is not active (multi-aggregated or fixed), it cannot be in any clique */
               if( cliquepartition[varseq[j]] == -1 && SCIPvarIsActive(tmpvars[varseq[j]]) )
               {
                  int k;

                  /* check if every variable in the actual clique is in clique with the new variable */
                  for( k = ncliquevars - 1; k >= 0; --k )
                  {
                     if( !SCIPvarsHaveCommonClique(tmpvars[varseq[j]], tmpvalues[varseq[j]], cliquevars[k],
                           cliquevalues[k], TRUE) )
                        break;
                  }

                  if( k == -1 )
                  {
                     /* put the variable into the same clique */
                     cliquepartition[varseq[j]] = cliquepartition[varseq[i]];
                     cliquevars[ncliquevars] = tmpvars[varseq[j]];
                     cliquevalues[ncliquevars] = tmpvalues[varseq[j]];
                     ++ncliquevars;
                  }
               }
            }
         }

         /* this clique is finished */
         ++(*ncliques);
      }
      assert(cliquepartition[varseq[i]] >= 0 && cliquepartition[varseq[i]] < i + 1);

      /* break if we reached the maximal number of comparisons */
      if( i * nvars > maxncliquevarscomp )
         break;
   }
   /* if we had too many variables fill up the cliquepartition and put each variable in a separate clique */
   for( ; i < nvars; ++i )
   {
      if( cliquepartition[varseq[i]] == -1 )
      {
         cliquepartition[varseq[i]] = *ncliques;
         ++(*ncliques);
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &sortkeys);
   SCIPfreeBufferArray(scip, &varseq);
   SCIPfreeBufferArray(scip, &tmpvars);
   SCIPfreeBufferArray(scip, &tmpvalues);
   SCIPfreeBufferArray(scip, &cliquevalues);
   SCIPfreeBufferArray(scip, &cliquevars);

   return SCIP_OKAY;
}

/** constructs sophisticated partition of knapsack variables into non-overlapping GUBs; current partition uses trivial GUBs */
static
SCIP_RETCODE GUBsetGetCliquePartition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure */
   SCIP_VAR**            vars,               /**< variables in the knapsack constraint */
   SCIP_Real*            solvals             /**< solution values of all knapsack variables */
   )
{
   int* cliquepartition;
   int* gubfirstvar;
   int ncliques;
   int currentgubconsidx;
   int newgubconsidx;
   int cliqueidx;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(gubset != NULL);
   assert(vars != NULL);

   nvars = gubset->nvars;
   assert(nvars >= 0);

   /* allocate temporary memory for clique partition */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquepartition, nvars) );

   /* compute sophisticated clique partition */
   SCIP_CALL( GUBsetCalcCliquePartition(scip, vars, nvars, cliquepartition, &ncliques, solvals) );

   /* allocate temporary memory for GUB set data structure */
   SCIP_CALL( SCIPallocBufferArray(scip, &gubfirstvar, ncliques) );

   /* translate GUB partition into GUB set data structure */
   for( i = 0; i < ncliques; i++ )
   {
      /* initialize first variable for every GUB */
      gubfirstvar[i] = -1;
   }
   /* move every knapsack variable into GUB defined by clique partition */
   for( i = 0; i < nvars; i++ )
   {
      assert(cliquepartition[i] >= 0);

      cliqueidx = cliquepartition[i];
      currentgubconsidx = gubset->gubconssidx[i];
      assert(gubset->gubconss[currentgubconsidx]->ngubvars == 1 );

      /* variable is first element in GUB constraint defined by clique partition */
      if( gubfirstvar[cliqueidx] == -1 )
      {
         /* corresponding GUB constraint in GUB set data structure was already constructed (as initial trivial GUB);
          * note: no assert for gubconssidx, because it can changed due to deleting empty GUBs in GUBsetMoveVar()
          */
         assert(gubset->gubvarsidx[i] == 0);
         assert(gubset->gubconss[gubset->gubconssidx[i]]->gubvars[gubset->gubvarsidx[i]] == i);

         /* remember the first variable found for the current GUB */
         gubfirstvar[cliqueidx] = i;
      }
      /* variable is additional element of GUB constraint defined by clique partition */
      else
      {
         assert(gubfirstvar[cliqueidx] >= 0 && gubfirstvar[cliqueidx] < i);

         /* move variable to GUB constraint defined by clique partition; index of this GUB constraint is given by the
          * first variable of this GUB constraint
          */
         newgubconsidx = gubset->gubconssidx[gubfirstvar[cliqueidx]];
         assert(newgubconsidx != currentgubconsidx); /* because initially every variable is in a different GUB */
         SCIP_CALL( GUBsetMoveVar(scip, gubset, vars, i, currentgubconsidx, newgubconsidx) );

         assert(gubset->gubconss[gubset->gubconssidx[i]]->gubvars[gubset->gubvarsidx[i]] == i);
      }
   }

#ifdef SCIP_DEBUG
   /* prints GUB set data structure */
   GUBsetPrint(scip, gubset, vars, solvals);
#endif

#ifndef NDEBUG
   /* checks consistency of GUB set data structure */
   SCIP_CALL( GUBsetCheck(scip, gubset, vars) );
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &gubfirstvar);
   SCIPfreeBufferArray(scip, &cliquepartition);

   return SCIP_OKAY;
}

/** gets a most violated cover C (\f$\sum_{j \in C} a_j > a_0\f$) for a given knapsack constraint \f$\sum_{j \in N} a_j x_j \leq a_0\f$
 *  taking into consideration the following fixing: \f$j \in C\f$, if \f$j \in N_1 = \{j \in N : x^*_j = 1\}\f$ and
 *  \f$j \in N \setminus C\f$, if \f$j \in N_0 = \{j \in N : x^*_j = 0\}\f$, if one exists.
 */
static
SCIP_RETCODE getCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< pointer to store cover variables */
   int*                  noncovervars,       /**< pointer to store noncover variables */
   int*                  ncovervars,         /**< pointer to store number of cover variables */
   int*                  nnoncovervars,      /**< pointer to store number of noncover variables */
   SCIP_Longint*         coverweight,        /**< pointer to store weight of cover */
   SCIP_Bool*            found,              /**< pointer to store whether a cover was found */
   SCIP_Bool             modtransused,       /**< should modified transformed separation problem be used to find cover */
   int*                  ntightened,         /**< pointer to store number of variables with tightened upper bound */
   SCIP_Bool*            fractional          /**< pointer to store whether the LP sol for knapsack vars is fractional */
   )
{
   SCIP_Longint* transweights;
   SCIP_Real* transprofits;
   SCIP_Longint transcapacity;
   SCIP_Longint fixedonesweight;
   SCIP_Longint itemsweight;
   SCIP_Bool infeasible;
   int* fixedones;
   int* fixedzeros;
   int* items;
   int nfixedones;
   int nfixedzeros;
   int nitems;
   int j;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars != NULL);
   assert(nnoncovervars != NULL);
   assert(coverweight != NULL);
   assert(found != NULL);
   assert(ntightened != NULL);
   assert(fractional != NULL);

   SCIPdebugMsg(scip, "   get cover for knapsack constraint\n");

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &transweights, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofits, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedones, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedzeros, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &items, nvars) );

   *found = FALSE;
   *ncovervars = 0;
   *nnoncovervars = 0;
   *coverweight = 0;
   *fractional = TRUE;

   /* gets the following sets
    *  N_1 = {j in N : x*_j = 1} (fixedones),
    *  N_0 = {j in N : x*_j = 0} (fixedzeros) and
    *  N\(N_0 & N_1) (items),
    * where x*_j is the solution value of variable x_j
    */
   nfixedones = 0;
   nfixedzeros = 0;
   nitems = 0;
   fixedonesweight = 0;
   itemsweight = 0;
   *ntightened = 0;
   for( j = 0; j < nvars; j++ )
   {
      assert(SCIPvarIsBinary(vars[j]));

      /* tightens upper bound of x_j if weight of x_j is greater than capacity of knapsack */
      if( weights[j] > capacity )
      {
         SCIP_CALL( SCIPtightenVarUb(scip, vars[j], 0.0, FALSE, &infeasible, NULL) );
         assert(!infeasible);
         (*ntightened)++;
         continue;
      }

      /* variable x_j has solution value one */
      if( SCIPisFeasEQ(scip, solvals[j], 1.0) )
      {
         fixedones[nfixedones] = j;
         nfixedones++;
         fixedonesweight += weights[j];
      }
      /* variable x_j has solution value zero */
      else if( SCIPisFeasEQ(scip, solvals[j], 0.0) )
      {
         fixedzeros[nfixedzeros] = j;
         nfixedzeros++;
      }
      /* variable x_j has fractional solution value */
      else
      {
         assert( SCIPisFeasGT(scip, solvals[j], 0.0) && SCIPisFeasLT(scip, solvals[j], 1.0) );
         items[nitems] = j;
         nitems++;
         itemsweight += weights[j];
      }
   }
   assert(nfixedones + nfixedzeros + nitems == nvars - (*ntightened));

   /* sets whether the LP solution x* for the knapsack variables is fractional; if it is not fractional we stop
    * the separation routine
    */
   assert(nitems >= 0);
   if( nitems == 0 )
   {
      *fractional = FALSE;
      goto TERMINATE;
   }
   assert(*fractional);

   /* transforms the traditional separation problem (under consideration of the following fixing:
    * z_j = 1 for all j in N_1, z_j = 0 for all j in N_0)
    *
    *   min sum_{j in N\(N_0 & N_1)} (1 - x*_j) z_j
    *       sum_{j in N\(N_0 & N_1)} a_j z_j >= (a_0 + 1) - sum_{j in N_1} a_j
    *                                    z_j in {0,1}, j in N\(N_0 & N_1)
    *
    * to a knapsack problem in maximization form by complementing the variables
    *
    * sum_{j in N\(N_0 & N_1)} (1 - x*_j) -
    *   max sum_{j in N\(N_0 & N_1)} (1 - x*_j) z_j
    *       sum_{j in N\(N_0 & N_1)} a_j z_j <= sum_{j in N\N_0} a_j - (a_0 + 1)
    *                                    z_j in {0,1}, j in N\(N_0 & N_1)
    */

   /* gets weight and profit of variables in transformed knapsack problem */
   for( j = 0; j < nitems; j++ )
   {
      transweights[j] = weights[items[j]];
      transprofits[j] = 1.0 - solvals[items[j]];
   }
   /* gets capacity of transformed knapsack problem */
   transcapacity = fixedonesweight + itemsweight - capacity - 1;

   /* if capacity of transformed knapsack problem is less than zero, there is no cover
    * (when variables fixed to zero are not used)
    */
   if( transcapacity < 0 )
   {
      assert(!(*found));
      goto TERMINATE;
   }

   if( modtransused )
   {
      /* transforms the modified separation problem (under consideration of the following fixing:
       * z_j = 1 for all j in N_1, z_j = 0 for all j in N_0)
       *
       *   min sum_{j in N\(N_0 & N_1)} (1 - x*_j) a_j z_j
       *       sum_{j in N\(N_0 & N_1)} a_j z_j >= (a_0 + 1) - sum_{j in N_1} a_j
       *                                    z_j in {0,1}, j in N\(N_0 & N_1)
       *
       * to a knapsack problem in maximization form by complementing the variables
       *
       * sum_{j in N\(N_0 & N_1)} (1 - x*_j) a_j -
       *   max sum_{j in N\(N_0 & N_1)} (1 - x*_j) a_j z_j
       *       sum_{j in N\(N_0 & N_1)} a_j z_j <= sum_{j in N\N_0} a_j - (a_0 + 1)
       *                                    z_j in {0,1}, j in N\(N_0 & N_1)
       */

      /* gets weight and profit of variables in modified transformed knapsack problem */
      for( j = 0; j < nitems; j++ )
      {
         transprofits[j] *= weights[items[j]];
         assert(SCIPisFeasPositive(scip, transprofits[j]));
      }
   }

   /* solves (modified) transformed knapsack problem approximately by solving the LP-relaxation of the (modified)
    * transformed knapsack problem using Dantzig's method and rounding down the solution.
    * let z* be the solution, then
    *   j in C,          if z*_j = 0 and
    *   i in N\C,        if z*_j = 1.
    */
   SCIP_CALL( SCIPsolveKnapsackApproximately(scip, nitems, transweights, transprofits, transcapacity, items,
         noncovervars, covervars, nnoncovervars, ncovervars, NULL) );
   /*assert(checkSolveKnapsack(scip, nitems, transweights, transprofits, items, weights, solvals, modtransused));*/

   /* constructs cover C (sum_{j in C} a_j > a_0) */
   for( j = 0; j < *ncovervars; j++ )
   {
      (*coverweight) += weights[covervars[j]];
   }

   /* adds all variables from N_1 to C */
   for( j = 0; j < nfixedones; j++ )
   {
      covervars[*ncovervars] = fixedones[j];
      (*ncovervars)++;
      (*coverweight) += weights[fixedones[j]];
   }

   /* adds all variables from N_0 to N\C */
   for( j = 0; j < nfixedzeros; j++ )
   {
      noncovervars[*nnoncovervars] = fixedzeros[j];
      (*nnoncovervars)++;
   }
   assert((*ncovervars) + (*nnoncovervars) == nvars - (*ntightened));
   assert((*coverweight) > capacity);
   *found = TRUE;

 TERMINATE:
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &items);
   SCIPfreeBufferArray(scip, &fixedzeros);
   SCIPfreeBufferArray(scip, &fixedones);
   SCIPfreeBufferArray(scip, &transprofits);
   SCIPfreeBufferArray(scip, &transweights);

   SCIPdebugMsg(scip, "   get cover for knapsack constraint -- end\n");

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks if minweightidx is set correctly
 */
static
SCIP_Bool checkMinweightidx(
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  covervars,          /**< pointer to store cover variables */
   int                   ncovervars,         /**< pointer to store number of cover variables */
   SCIP_Longint          coverweight,        /**< pointer to store weight of cover */
   int                   minweightidx,       /**< index of variable in cover variables with minimum weight */
   int                   j                   /**< current index in cover variables */
   )
{
   SCIP_Longint minweight;
   int i;

   assert(weights != NULL);
   assert(covervars != NULL);
   assert(ncovervars > 0);

   minweight = weights[covervars[minweightidx]];

   /* checks if all cover variables before index j have weight greater than minweight */
   for( i = 0; i < j; i++ )
   {
      assert(weights[covervars[i]] > minweight);
      if( weights[covervars[i]] <= minweight )
         return FALSE;
   }

   /* checks if all variables before index j cannot be removed, i.e. i cannot be the next minweightidx */
   for( i = 0; i < j; i++ )
   {
      assert(coverweight - weights[covervars[i]] <= capacity);
      if( coverweight - weights[covervars[i]] > capacity )
         return FALSE;
   }
   return TRUE;
}
#endif


/** gets partition \f$(C_1,C_2)\f$ of minimal cover \f$C\f$, i.e. \f$C_1 \cup C_2 = C\f$ and \f$C_1 \cap C_2 = \emptyset\f$,
 *  with \f$C_1\f$ not empty; chooses partition as follows \f$C_2 = \{ j \in C : x^*_j = 1 \}\f$ and \f$C_1 = C \setminus C_2\f$
 */
static
void getPartitionCovervars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables */
   int                   ncovervars,         /**< number of cover variables */
   int*                  varsC1,             /**< pointer to store variables in C1 */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  nvarsC1,            /**< pointer to store number of variables in C1 */
   int*                  nvarsC2             /**< pointer to store number of variables in C2 */
   )
{
   int j;

   assert(scip != NULL);
   assert(ncovervars >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(varsC1 != NULL);
   assert(varsC2 != NULL);
   assert(nvarsC1 != NULL);
   assert(nvarsC2 != NULL);

   *nvarsC1 = 0;
   *nvarsC2 = 0;
   for( j = 0; j < ncovervars; j++ )
   {
      assert(SCIPisFeasGT(scip, solvals[covervars[j]], 0.0));

      /* variable has solution value one */
      if( SCIPisGE(scip, solvals[covervars[j]], 1.0) )
      {
         varsC2[*nvarsC2] = covervars[j];
         (*nvarsC2)++;
      }
      /* variable has solution value less than one */
      else
      {
         assert(SCIPisLT(scip, solvals[covervars[j]], 1.0));
         varsC1[*nvarsC1] = covervars[j];
         (*nvarsC1)++;
      }
   }
   assert((*nvarsC1) + (*nvarsC2) == ncovervars);
}

/** changes given partition (C_1,C_2) of minimal cover C, if |C1| = 1, by moving one and two (if possible) variables from
 *  C2 to C1 if |C1| = 1 and |C1| = 0, respectively.
 */
static
SCIP_RETCODE changePartitionCovervars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   int*                  varsC1,             /**< pointer to store variables in C1 */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  nvarsC1,            /**< pointer to store number of variables in C1 */
   int*                  nvarsC2             /**< pointer to store number of variables in C2 */
   )
{
   SCIP_Real* sortkeysC2;
   int j;

   assert(*nvarsC1 >= 0 && *nvarsC1 <= 1);
   assert(*nvarsC2 > 0);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC2, *nvarsC2) );

   /* sorts variables in C2 such that a_1 >= .... >= a_|C2| */
   for( j = 0; j < *nvarsC2; j++ )
      sortkeysC2[j] = (SCIP_Real) weights[varsC2[j]];
   SCIPsortDownRealInt(sortkeysC2, varsC2, *nvarsC2);

   /* adds one or two variable from C2 with smallest weight to C1 and removes them from C2 */
   assert(*nvarsC2 == 1 || weights[varsC2[(*nvarsC2)-1]] <= weights[varsC2[(*nvarsC2)-2]]);
   while( *nvarsC1 < 2 && *nvarsC2 > 0 )
   {
      varsC1[*nvarsC1] = varsC2[(*nvarsC2)-1];
      (*nvarsC1)++;
      (*nvarsC2)--;
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeysC2);

   return SCIP_OKAY;
}

/** changes given partition (C_1,C_2) of feasible set C, if |C1| = 1, by moving one variable from C2 to C1 */
static
SCIP_RETCODE changePartitionFeasiblesetvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   int*                  varsC1,             /**< pointer to store variables in C1 */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  nvarsC1,            /**< pointer to store number of variables in C1 */
   int*                  nvarsC2             /**< pointer to store number of variables in C2 */
   )
{
   SCIP_Real* sortkeysC2;
   int j;

   assert(*nvarsC1 >= 0 && *nvarsC1 <= 1);
   assert(*nvarsC2 > 0);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC2, *nvarsC2) );

   /* sorts variables in C2 such that a_1 >= .... >= a_|C2| */
   for( j = 0; j < *nvarsC2; j++ )
      sortkeysC2[j] = (SCIP_Real) weights[varsC2[j]];
   SCIPsortDownRealInt(sortkeysC2, varsC2, *nvarsC2);

   /* adds variable from C2 with smallest weight to C1 and removes it from C2 */
   assert(*nvarsC2 == 1 || weights[varsC2[(*nvarsC2)-1]] <= weights[varsC2[(*nvarsC2)-2]]);
   varsC1[*nvarsC1] = varsC2[(*nvarsC2)-1];
   (*nvarsC1)++;
   (*nvarsC2)--;

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeysC2);

   return SCIP_OKAY;
}


/** gets partition \f$(F,R)\f$ of \f$N \setminus C\f$ where \f$C\f$ is a minimal cover, i.e. \f$F \cup R = N \setminus C\f$
 *  and \f$F \cap R = \emptyset\f$; chooses partition as follows \f$R = \{ j \in N \setminus C : x^*_j = 0 \}\f$ and
 *  \f$F = (N \setminus C) \setminus F\f$
 */
static
void getPartitionNoncovervars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  noncovervars,       /**< noncover variables */
   int                   nnoncovervars,      /**< number of noncover variables */
   int*                  varsF,              /**< pointer to store variables in F */
   int*                  varsR,              /**< pointer to store variables in R */
   int*                  nvarsF,             /**< pointer to store number of variables in F */
   int*                  nvarsR              /**< pointer to store number of variables in R */
   )
{
   int j;

   assert(scip != NULL);
   assert(nnoncovervars >= 0);
   assert(solvals != NULL);
   assert(noncovervars != NULL);
   assert(varsF != NULL);
   assert(varsR != NULL);
   assert(nvarsF != NULL);
   assert(nvarsR != NULL);

   *nvarsF = 0;
   *nvarsR = 0;

   for( j = 0; j < nnoncovervars; j++ )
   {
      /* variable has solution value zero */
      if( SCIPisFeasEQ(scip, solvals[noncovervars[j]], 0.0) )
      {
         varsR[*nvarsR] = noncovervars[j];
         (*nvarsR)++;
      }
      /* variable has solution value greater than zero */
      else
      {
         assert(SCIPisFeasGT(scip, solvals[noncovervars[j]], 0.0));
         varsF[*nvarsF] = noncovervars[j];
         (*nvarsF)++;
      }
   }
   assert((*nvarsF) + (*nvarsR) == nnoncovervars);
}

/** sorts variables in F, C_2, and R according to the second level lifting sequence that will be used in the sequential
 *  lifting procedure
 */
static
SCIP_RETCODE getLiftingSequence(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   int*                  varsF,              /**< pointer to store variables in F */
   int*                  varsC2,             /**< pointer to store variables in C2 */
   int*                  varsR,              /**< pointer to store variables in R */
   int                   nvarsF,             /**< number of variables in F */
   int                   nvarsC2,            /**< number of variables in C2 */
   int                   nvarsR              /**< number of variables in R */
   )
{
   SORTKEYPAIR** sortkeypairsF;
   SORTKEYPAIR* sortkeypairsFstore;
   SCIP_Real* sortkeysC2;
   SCIP_Real* sortkeysR;
   int j;

   assert(scip != NULL);
   assert(solvals != NULL);
   assert(weights != NULL);
   assert(varsF != NULL);
   assert(varsC2 != NULL);
   assert(varsR != NULL);
   assert(nvarsF >= 0);
   assert(nvarsC2 >= 0);
   assert(nvarsR >= 0);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeypairsF, nvarsF) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeypairsFstore, nvarsF) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC2, nvarsC2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysR, nvarsR) );

   /* gets sorting key for variables in F corresponding to the following lifting sequence
    *  sequence 1: non-increasing absolute difference between x*_j and the value the variable is fixed to, i.e.
    *              x*_1 >= x*_2 >= ... >= x*_|F|
    * in case of equality uses
    *  sequence 4: non-increasing a_j, i.e. a_1 >= a_2 >= ... >= a_|C_2|
    */
   for( j = 0; j < nvarsF; j++ )
   {
      sortkeypairsF[j] = &(sortkeypairsFstore[j]);
      sortkeypairsF[j]->key1 = solvals[varsF[j]];
      sortkeypairsF[j]->key2 = (SCIP_Real) weights[varsF[j]];
   }

   /* gets sorting key for variables in C_2 corresponding to the following lifting sequence
    *  sequence 4: non-increasing a_j, i.e. a_1 >= a_2 >= ... >= a_|C_2|
    */
   for( j = 0; j < nvarsC2; j++ )
      sortkeysC2[j] = (SCIP_Real) weights[varsC2[j]];

   /* gets sorting key for variables in R corresponding to the following lifting sequence
    *  sequence 4: non-increasing a_j, i.e. a_1 >= a_2 >= ... >= a_|R|
    */
   for( j = 0; j < nvarsR; j++ )
      sortkeysR[j] = (SCIP_Real) weights[varsR[j]];

   /* sorts F, C2 and R */
   if( nvarsF > 0 )
   {
      SCIPsortDownPtrInt((void**)sortkeypairsF, varsF, compSortkeypairs, nvarsF);
   }
   if( nvarsC2 > 0 )
   {
      SCIPsortDownRealInt(sortkeysC2, varsC2, nvarsC2);
   }
   if( nvarsR > 0)
   {
      SCIPsortDownRealInt(sortkeysR, varsR, nvarsR);
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeysR);
   SCIPfreeBufferArray(scip, &sortkeysC2);
   SCIPfreeBufferArray(scip, &sortkeypairsFstore);
   SCIPfreeBufferArray(scip, &sortkeypairsF);

   return SCIP_OKAY;
}

/** categorizes GUBs of knapsack GUB partion into GOC1, GNC1, GF, GC2, and GR and computes a lifting sequence of the GUBs
 *  for the sequential GUB wise lifting procedure
 */
static
SCIP_RETCODE getLiftingSequenceGUB(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure */
   SCIP_Real*            solvals,            /**< solution values of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   int*                  varsC1,             /**< variables in C1 */
   int*                  varsC2,             /**< variables in C2 */
   int*                  varsF,              /**< variables in F */
   int*                  varsR,              /**< variables in R */
   int                   nvarsC1,            /**< number of variables in C1 */
   int                   nvarsC2,            /**< number of variables in C2 */
   int                   nvarsF,             /**< number of variables in F */
   int                   nvarsR,             /**< number of variables in R */
   int*                  gubconsGC1,         /**< pointer to store GUBs in GC1(GNC1+GOC1) */
   int*                  gubconsGC2,         /**< pointer to store GUBs in GC2 */
   int*                  gubconsGFC1,        /**< pointer to store GUBs in GFC1(GNC1+GF) */
   int*                  gubconsGR,          /**< pointer to store GUBs in GR */
   int*                  ngubconsGC1,        /**< pointer to store number of GUBs in GC1(GNC1+GOC1) */
   int*                  ngubconsGC2,        /**< pointer to store number of GUBs in GC2 */
   int*                  ngubconsGFC1,       /**< pointer to store number of GUBs in GFC1(GNC1+GF) */
   int*                  ngubconsGR,         /**< pointer to store number of GUBs in GR */
   int*                  ngubconscapexceed,  /**< pointer to store number of GUBs with only capacity exceeding variables */
   int*                  maxgubvarssize      /**< pointer to store the maximal size of GUB constraints */
   )
{
#if 0 /* not required */
   SORTKEYPAIR** sortkeypairsF;
#endif
   SORTKEYPAIR** sortkeypairsGFC1;
   SORTKEYPAIR* sortkeypairsGFC1store;
   SCIP_Real* sortkeysC1;
   SCIP_Real* sortkeysC2;
   SCIP_Real* sortkeysR;
   int* nC1varsingubcons;
   int var;
   int gubconsidx;
   int varidx;
   int ngubconss;
   int ngubconsGOC1;
   int targetvar;
   int nvarsprocessed;
   int i;
   int j;

#if GUBSPLITGNC1GUBS
   SCIP_Bool gubconswithF;
   int origngubconss;
   origngubconss = gubset->ngubconss;
#endif

   assert(scip != NULL);
   assert(gubset != NULL);
   assert(solvals != NULL);
   assert(weights != NULL);
   assert(varsC1 != NULL);
   assert(varsC2 != NULL);
   assert(varsF != NULL);
   assert(varsR != NULL);
   assert(nvarsC1 > 0);
   assert(nvarsC2 >= 0);
   assert(nvarsF >= 0);
   assert(nvarsR >= 0);
   assert(gubconsGC1 != NULL);
   assert(gubconsGC2 != NULL);
   assert(gubconsGFC1 != NULL);
   assert(gubconsGR != NULL);
   assert(ngubconsGC1 != NULL);
   assert(ngubconsGC2 != NULL);
   assert(ngubconsGFC1 != NULL);
   assert(ngubconsGR != NULL);
   assert(maxgubvarssize != NULL);

   ngubconss = gubset->ngubconss;
   nvarsprocessed = 0;
   ngubconsGOC1 = 0;

   /* GUBs are categorized into different types according to the variables in volved
    * - GOC1:  involves variables in C1 only           -- no C2, R, F
    * - GNC1:  involves variables in C1 and F (and R)  -- no C2
    * - GF:    involves variables in F  (and R) only   -- no C1, C2
    * - GC2:   involves variables in C2 only           -- no C1, R, F
    * - GR:    involves variables in R  only           -- no C1, C2, F
    * which requires splitting GUBs in case they include variable in F and R.
    *
    * afterwards all GUBs (except GOC1 GUBs, which we do not need to lift) are sorted by a two level lifting sequence.
    * - first  ordering level is: GFC1 (GNC1+GF), GC2, and GR.
    * - second ordering level is
    *    GFC1:   non-increasing number of variables in F and non-increasing max{x*_k : k in GFC1_j} in case of equality
    *    GC2:    non-increasing max{ a_k : k in GC2_j}; note that |GFC2_j| = 1
    *    GR:     non-increasing max{ a_k : k in GR_j}
    *
    * in additon, another GUB union, which is helpful for the lifting procedure, is formed
    * - GC1:   GUBs of category GOC1 and GNC1
    * with second ordering level non-decreasing min{ a_k : k in GC1_j };
    * note that min{ a_k : k in GC1_j } always comes from the first variable in the GUB
    */

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC1, nvarsC1) );
#if 0 /* not required */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeypairsF, nvarsF) );
#endif
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysC2, nvarsC2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeysR, nvarsR) );


   /* to get the GUB lifting sequence, we first sort all variables in F, C2, and R
    * - F:      non-increasing x*_j and non-increasing a_j in case of equality
    * - C2:     non-increasing a_j
    * - R:      non-increasing a_j
    * furthermore, sort C1 variables as needed for initializing the minweight table (non-increasing a_j).
    */

   /* gets sorting key for variables in C1 corresponding to the following ordering
    *  non-decreasing a_j, i.e. a_1 <= a_2 <= ... <= a_|C_1|
    */
   for( j = 0; j < nvarsC1; j++ )
   {
      /* gets sortkeys */
      sortkeysC1[j] = (SCIP_Real) weights[varsC1[j]];

      /* update status of variable in its gub constraint */
      gubconsidx = gubset->gubconssidx[varsC1[j]];
      varidx = gubset->gubvarsidx[varsC1[j]];
      gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] = GUBVARSTATUS_BELONGSTOSET_C1;
   }

   /* gets sorting key for variables in F corresponding to the following ordering
    *  non-increasing x*_j, i.e., x*_1 >= x*_2 >= ... >= x*_|F|, and
    *  non-increasing a_j,  i.e., a_1  >= a_2  >= ... >= a_|F| in case of equality
    * and updates status of each variable in F in GUB set data structure
    */
   for( j = 0; j < nvarsF; j++ )
   {
#if 0 /* not required */
      /* gets sortkeys */
      SCIP_CALL( SCIPallocBuffer(scip, &sortkeypairsF[j]) );
      sortkeypairsF[j]->key1 = solvals[varsF[j]];
      sortkeypairsF[j]->key2 = (SCIP_Real) weights[varsF[j]];
#endif

      /* update status of variable in its gub constraint */
      gubconsidx = gubset->gubconssidx[varsF[j]];
      varidx = gubset->gubvarsidx[varsF[j]];
      gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] = GUBVARSTATUS_BELONGSTOSET_F;
   }

   /* gets sorting key for variables in C2 corresponding to the following ordering
    *  non-increasing a_j,  i.e., a_1  >= a_2  >= ... >= a_|C2|
    * and updates status of each variable in F in GUB set data structure
    */
   for( j = 0; j < nvarsC2; j++ )
   {
      /* gets sortkeys */
      sortkeysC2[j] = (SCIP_Real) weights[varsC2[j]];

      /* update status of variable in its gub constraint */
      gubconsidx = gubset->gubconssidx[varsC2[j]];
      varidx = gubset->gubvarsidx[varsC2[j]];
      gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] = GUBVARSTATUS_BELONGSTOSET_C2;
   }

   /* gets sorting key for variables in R corresponding to the following ordering
    *  non-increasing a_j,  i.e., a_1  >= a_2  >= ... >= a_|R|
    * and updates status of each variable in F in GUB set data structure
    */
   for( j = 0; j < nvarsR; j++ )
   {
      /* gets sortkeys */
      sortkeysR[j] = (SCIP_Real) weights[varsR[j]];

      /* update status of variable in its gub constraint */
      gubconsidx = gubset->gubconssidx[varsR[j]];
      varidx = gubset->gubvarsidx[varsR[j]];
      gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] = GUBVARSTATUS_BELONGSTOSET_R;
   }

   /* sorts C1, F, C2 and R */
   if( nvarsC1 > 0 )
   {
      SCIPsortRealInt(sortkeysC1, varsC1, nvarsC1);
   }
#if 0 /* not required */
   if( nvarsF > 0 )
   {
      SCIPsortDownPtrInt((void**)sortkeypairsF, varsF, compSortkeypairs, nvarsF);
   }
#endif
   if( nvarsC2 > 0 )
   {
      SCIPsortDownRealInt(sortkeysC2, varsC2, nvarsC2);
   }
   if( nvarsR > 0)
   {
      SCIPsortDownRealInt(sortkeysR, varsR, nvarsR);
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeysR);
   SCIPfreeBufferArray(scip, &sortkeysC2);
#if 0 /* not required */
   for( j = nvarsF-1; j >= 0; j-- )
      SCIPfreeBuffer(scip, &sortkeypairsF[j]);
   SCIPfreeBufferArray(scip, &sortkeypairsF);
#endif
   SCIPfreeBufferArray(scip, &sortkeysC1);

   /* allocate and initialize temporary memory for sorting GUB constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeypairsGFC1, ngubconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeypairsGFC1store, ngubconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nC1varsingubcons, ngubconss) );
   BMSclearMemoryArray(nC1varsingubcons, ngubconss);
   for( i = 0; i < ngubconss; i++)
   {
      sortkeypairsGFC1[i] = &(sortkeypairsGFC1store[i]);
      sortkeypairsGFC1[i]->key1 = 0.0;
      sortkeypairsGFC1[i]->key2 = 0.0;
   }
   *ngubconsGC1 = 0;
   *ngubconsGC2 = 0;
   *ngubconsGFC1 = 0;
   *ngubconsGR = 0;
   *ngubconscapexceed = 0;
   *maxgubvarssize = 0;

#ifndef NDEBUG
   for( i = 0; i < gubset->ngubconss; i++ )
      assert(gubset->gubconsstatus[i] == GUBCONSSTATUS_UNINITIAL);
#endif

   /* stores GUBs of group GC1 (GOC1+GNC1) and part of the GUBs of group GFC1 (GNC1 GUBs) and sorts variables in these GUBs
    * s.t. C1 variables come first (will automatically be sorted by non-decreasing weight).
    * gets sorting keys for GUBs of type GFC1 corresponding to the following ordering
    *    non-increasing number of variables in F, and
    *    non-increasing max{x*_k : k in GFC1_j} in case of equality
    */
   for( i = 0; i < nvarsC1; i++ )
   {
      int nvarsC1capexceed;

      nvarsC1capexceed = 0;

      var = varsC1[i];
      gubconsidx = gubset->gubconssidx[var];
      varidx = gubset->gubvarsidx[var];

      assert(gubconsidx >= 0 && gubconsidx < ngubconss);
      assert(gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] == GUBVARSTATUS_BELONGSTOSET_C1);

      /* current C1 variable is put to the front of its GUB where C1 part is stored by non-decreasing weigth;
       * note that variables in C1 are already sorted by non-decreasing weigth
       */
      targetvar = gubset->gubconss[gubconsidx]->gubvars[nC1varsingubcons[gubconsidx]];
      GUBsetSwapVars(scip, gubset, var, targetvar);
      nC1varsingubcons[gubconsidx]++;

      /* the GUB was already handled (status set and stored in its group) by another variable of the GUB */
      if( gubset->gubconsstatus[gubconsidx] != GUBCONSSTATUS_UNINITIAL )
      {
         assert(gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GOC1
            || gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1);
         continue;
      }

      /* determine the status of the current GUB constraint, GOC1 or GNC1; GUBs involving R variables are split into
       * GOC1/GNC1 and GF, if wanted. also update sorting key if GUB is of type GFC1 (GNC1)
       */
#if GUBSPLITGNC1GUBS
      gubconswithF = FALSE;
#endif
      for( j = 0; j < gubset->gubconss[gubconsidx]->ngubvars; j++ )
      {
         assert(gubset->gubconss[gubconsidx]->gubvarsstatus[j] != GUBVARSTATUS_BELONGSTOSET_C2);

         /* C1-variable: update number of C1/capacity exceeding variables */
         if( gubset->gubconss[gubconsidx]->gubvarsstatus[j] == GUBVARSTATUS_BELONGSTOSET_C1 )
         {
            nvarsC1capexceed++;
            nvarsprocessed++;
         }
         /* F-variable: update sort key (number of F variables in GUB) of corresponding GFC1-GUB */
         else if( gubset->gubconss[gubconsidx]->gubvarsstatus[j] == GUBVARSTATUS_BELONGSTOSET_F )
         {
#if GUBSPLITGNC1GUBS
	    gubconswithF = TRUE;
#endif
	    sortkeypairsGFC1[*ngubconsGFC1]->key1 += 1.0;

            if( solvals[gubset->gubconss[gubconsidx]->gubvars[j]] > sortkeypairsGFC1[*ngubconsGFC1]->key2 )
               sortkeypairsGFC1[*ngubconsGFC1]->key2 = solvals[gubset->gubconss[gubconsidx]->gubvars[j]];
         }
         else if( gubset->gubconss[gubconsidx]->gubvarsstatus[j] == GUBVARSTATUS_CAPACITYEXCEEDED )
         {
            nvarsC1capexceed++;
         }
         else
            assert(gubset->gubconss[gubconsidx]->gubvarsstatus[j] == GUBVARSTATUS_BELONGSTOSET_R);
      }

      /* update set of GC1 GUBs */
      gubconsGC1[*ngubconsGC1] = gubconsidx;
      (*ngubconsGC1)++;

      /* update maximum size of all GUB constraints */
      if( gubset->gubconss[gubconsidx]->gubvarssize > *maxgubvarssize )
	 *maxgubvarssize = gubset->gubconss[gubconsidx]->gubvarssize;

      /* set status of GC1-GUB (GOC1 or GNC1) and update set of GFC1 GUBs */
      if( nvarsC1capexceed == gubset->gubconss[gubconsidx]->ngubvars )
      {
         gubset->gubconsstatus[gubconsidx] = GUBCONSSTATUS_BELONGSTOSET_GOC1;
	 ngubconsGOC1++;
      }
      else
      {
#if GUBSPLITGNC1GUBS
         /* only variables in C1 and R -- no in F: GUB will be split into GR and GOC1 GUBs */
	 if( !gubconswithF )
	 {
	    GUBVARSTATUS movevarstatus;

	    assert(gubset->ngubconss < gubset->nvars);

            /* create a new GUB for GR part of splitting */
	    SCIP_CALL( GUBconsCreate(scip, &gubset->gubconss[gubset->ngubconss]) );
	    gubset->ngubconss++;
	    ngubconss = gubset->ngubconss;

            /* fill GR with R variables in current GUB */
	    for( j = gubset->gubconss[gubconsidx]->ngubvars-1; j >= 0; j-- )
	    {
	        movevarstatus = gubset->gubconss[gubconsidx]->gubvarsstatus[j];
		if( movevarstatus != GUBVARSTATUS_BELONGSTOSET_C1 )
		{
		   assert(movevarstatus == GUBVARSTATUS_BELONGSTOSET_R || movevarstatus == GUBVARSTATUS_CAPACITYEXCEEDED);
		   SCIP_CALL( GUBsetMoveVar(scip, gubset, vars, gubset->gubconss[gubconsidx]->gubvars[j],
                         gubconsidx, ngubconss-1) );
		   gubset->gubconss[ngubconss-1]->gubvarsstatus[gubset->gubconss[ngubconss-1]->ngubvars-1] =
                      movevarstatus;
		}
	    }

	    gubset->gubconsstatus[gubconsidx] = GUBCONSSTATUS_BELONGSTOSET_GOC1;
	    ngubconsGOC1++;

	    gubset->gubconsstatus[ngubconss-1] = GUBCONSSTATUS_BELONGSTOSET_GR;
	    gubconsGR[*ngubconsGR] = ngubconss-1;
	    (*ngubconsGR)++;
	 }
         /* variables in C1, F, and maybe R: GNC1 GUB */
	 else
	 {
	    assert(gubconswithF);

	    gubset->gubconsstatus[gubconsidx] = GUBCONSSTATUS_BELONGSTOSET_GNC1;
	    gubconsGFC1[*ngubconsGFC1] = gubconsidx;
	    (*ngubconsGFC1)++;
	 }
#else
	 gubset->gubconsstatus[gubconsidx] = GUBCONSSTATUS_BELONGSTOSET_GNC1;
	 gubconsGFC1[*ngubconsGFC1] = gubconsidx;
	 (*ngubconsGFC1)++;
#endif
      }
   }

   /* stores GUBs of group GC2 (only trivial GUBs); sorting is not required because the C2 variables (which we loop over)
    * are already sorted correctly
    */
   for( i = 0; i < nvarsC2; i++ )
   {
      var = varsC2[i];
      gubconsidx = gubset->gubconssidx[var];
      varidx = gubset->gubvarsidx[var];

      assert(gubconsidx >= 0 && gubconsidx < ngubconss);
      assert(gubset->gubconss[gubconsidx]->ngubvars == 1);
      assert(varidx == 0);
      assert(gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] == GUBVARSTATUS_BELONGSTOSET_C2);
      assert(gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_UNINITIAL);

      /* set status of GC2 GUB */
      gubset->gubconsstatus[gubconsidx] = GUBCONSSTATUS_BELONGSTOSET_GC2;

      /* update group of GC2 GUBs */
      gubconsGC2[*ngubconsGC2] = gubconsidx;
      (*ngubconsGC2)++;

      /* update maximum size of all GUB constraints */
      if( gubset->gubconss[gubconsidx]->gubvarssize > *maxgubvarssize )
	 *maxgubvarssize = gubset->gubconss[gubconsidx]->gubvarssize;

      nvarsprocessed++;
   }

   /* stores remaining part of the GUBs of group GFC1 (GF GUBs) and gets GUB sorting keys corresp. to following ordering
    *    non-increasing number of variables in F, and
    *    non-increasing max{x*_k : k in GFC1_j} in case of equality
    */
   for( i = 0; i < nvarsF; i++ )
   {
      var = varsF[i];
      gubconsidx = gubset->gubconssidx[var];
      varidx = gubset->gubvarsidx[var];

      assert(gubconsidx >= 0 && gubconsidx < ngubconss);
      assert(gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] == GUBVARSTATUS_BELONGSTOSET_F);

      nvarsprocessed++;

      /* the GUB was already handled (status set and stored in its group) by another variable of the GUB */
      if( gubset->gubconsstatus[gubconsidx] != GUBCONSSTATUS_UNINITIAL )
      {
	 assert(gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GF
	      || gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1);
         continue;
      }

      /* set status of GF GUB */
      gubset->gubconsstatus[gubconsidx] = GUBCONSSTATUS_BELONGSTOSET_GF;

      /* update sorting key of corresponding GFC1 GUB */
      for( j = 0; j < gubset->gubconss[gubconsidx]->ngubvars; j++ )
      {
         assert(gubset->gubconss[gubconsidx]->gubvarsstatus[j] != GUBVARSTATUS_BELONGSTOSET_C2
            && gubset->gubconss[gubconsidx]->gubvarsstatus[j] != GUBVARSTATUS_BELONGSTOSET_C1);

         /* F-variable: update sort key (number of F variables in GUB) of corresponding GFC1-GUB */
         if( gubset->gubconss[gubconsidx]->gubvarsstatus[j] == GUBVARSTATUS_BELONGSTOSET_F )
         {
            sortkeypairsGFC1[*ngubconsGFC1]->key1 += 1.0;

            if( solvals[gubset->gubconss[gubconsidx]->gubvars[j]] > sortkeypairsGFC1[*ngubconsGFC1]->key2 )
               sortkeypairsGFC1[*ngubconsGFC1]->key2 = solvals[gubset->gubconss[gubconsidx]->gubvars[j]];
         }
      }

      /* update set of GFC1 GUBs */
      gubconsGFC1[*ngubconsGFC1] = gubconsidx;
      (*ngubconsGFC1)++;

      /* update maximum size of all GUB constraints */
      if( gubset->gubconss[gubconsidx]->gubvarssize > *maxgubvarssize )
         *maxgubvarssize = gubset->gubconss[gubconsidx]->gubvarssize;
   }

   /* stores GUBs of group GR; sorting is not required because the R variables (which we loop over) are already sorted
    * correctly
    */
   for( i = 0; i < nvarsR; i++ )
   {
      var = varsR[i];
      gubconsidx = gubset->gubconssidx[var];
      varidx = gubset->gubvarsidx[var];

      assert(gubconsidx >= 0 && gubconsidx < ngubconss);
      assert(gubset->gubconss[gubconsidx]->gubvarsstatus[varidx] == GUBVARSTATUS_BELONGSTOSET_R);

      nvarsprocessed++;

      /* the GUB was already handled (status set and stored in its group) by another variable of the GUB */
      if( gubset->gubconsstatus[gubconsidx] != GUBCONSSTATUS_UNINITIAL )
      {
	 assert(gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GR
	      || gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GF
	      || gubset->gubconsstatus[gubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1);
         continue;
      }

      /* set status of GR GUB */
      gubset->gubconsstatus[gubconsidx] = GUBCONSSTATUS_BELONGSTOSET_GR;

      /* update set of GR GUBs */
      gubconsGR[*ngubconsGR] = gubconsidx;
      (*ngubconsGR)++;

      /* update maximum size of all GUB constraints */
      if( gubset->gubconss[gubconsidx]->gubvarssize > *maxgubvarssize )
         *maxgubvarssize = gubset->gubconss[gubconsidx]->gubvarssize;
   }
   assert(nvarsprocessed == nvarsC1 + nvarsC2 + nvarsF + nvarsR);

   /* update number of GUBs with only capacity exceeding variables (will not be used for lifting) */
   (*ngubconscapexceed) =  ngubconss - (ngubconsGOC1 + (*ngubconsGC2) + (*ngubconsGFC1) + (*ngubconsGR));
   assert(*ngubconscapexceed >= 0);
#ifndef NDEBUG
   {
      int check;

      check = 0;

      /* remaining not handled GUBs should only contain capacity exceeding variables */
      for( i = 0; i < ngubconss; i++ )
      {
         if( gubset->gubconsstatus[i] ==  GUBCONSSTATUS_UNINITIAL )
            check++;
      }
      assert(check == *ngubconscapexceed);
   }
#endif

   /* sort GFCI GUBs according to computed sorting keys */
   if( (*ngubconsGFC1) > 0 )
   {
      SCIPsortDownPtrInt((void**)sortkeypairsGFC1, gubconsGFC1, compSortkeypairs, (*ngubconsGFC1));
   }

   /* free temporary memory */
#if GUBSPLITGNC1GUBS
   ngubconss = origngubconss;
#endif
   SCIPfreeBufferArray(scip, &nC1varsingubcons);
   SCIPfreeBufferArray(scip, &sortkeypairsGFC1store);
   SCIPfreeBufferArray(scip, &sortkeypairsGFC1);

   return SCIP_OKAY;
}

/** enlarges minweight table to at least the given length */
static
SCIP_RETCODE enlargeMinweights(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint**        minweightsptr,      /**< pointer to minweights table */
   int*                  minweightslen,      /**< pointer to store number of entries in minweights table (incl. z=0) */
   int*                  minweightssize,     /**< pointer to current size of minweights table */
   int                   newlen              /**< new length of minweights table */
   )
{
   int j;

   assert(minweightsptr != NULL);
   assert(*minweightsptr != NULL);
   assert(minweightslen != NULL);
   assert(*minweightslen >= 0);
   assert(minweightssize != NULL);
   assert(*minweightssize >= 0);

   if( newlen > *minweightssize )
   {
      int newsize;

      /* reallocate table memory */
      newsize = SCIPcalcMemGrowSize(scip, newlen);
      SCIP_CALL( SCIPreallocBufferArray(scip, minweightsptr, newsize) );
      *minweightssize = newsize;
   }
   assert(newlen <= *minweightssize);

   /* initialize new elements */
   for( j = *minweightslen; j < newlen; ++j )
      (*minweightsptr)[j] = SCIP_LONGINT_MAX;
   *minweightslen = newlen;

   return SCIP_OKAY;
}

/** lifts given inequality
 *    sum_{j in M_1} x_j <= alpha_0
 *  valid for
 *    S^0 = { x in {0,1}^|M_1| : sum_{j in M_1} a_j x_j <= a_0 - sum_{j in M_2} a_j }
 *  to a valid inequality
 *    sum_{j in M_1} x_j + sum_{j in F} alpha_j x_j + sum_{j in M_2} alpha_j x_j + sum_{j in R} alpha_j x_j
 *    <= alpha_0 + sum_{j in M_2} alpha_j
 *  for
 *    S = { x in {0,1}^|N| : sum_{j in N} a_j x_j <= a_0 };
 *  uses sequential up-lifting for the variables in F, sequential down-lifting for the variable in M_2, and
 *  sequential up-lifting for the variables in R; procedure can be used to strengthen minimal cover inequalities and
 *  extended weight inequalities.
 */
static
SCIP_RETCODE sequentialUpAndDownLifting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  varsM1,             /**< variables in M_1 */
   int*                  varsM2,             /**< variables in M_2 */
   int*                  varsF,              /**< variables in F */
   int*                  varsR,              /**< variables in R */
   int                   nvarsM1,            /**< number of variables in M_1 */
   int                   nvarsM2,            /**< number of variables in M_2 */
   int                   nvarsF,             /**< number of variables in F */
   int                   nvarsR,             /**< number of variables in R */
   int                   alpha0,             /**< rights hand side of given valid inequality */
   int*                  liftcoefs,          /**< pointer to store lifting coefficient of vars in knapsack constraint */
   SCIP_Real*            cutact,             /**< pointer to store activity of lifted valid inequality */
   int*                  liftrhs             /**< pointer to store right hand side of the lifted valid inequality */
   )
{
   SCIP_Longint* minweights;
   SCIP_Real* sortkeys;
   SCIP_Longint fixedonesweight;
   int minweightssize;
   int minweightslen;
   int j;
   int w;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(varsM1 != NULL);
   assert(varsM2 != NULL);
   assert(varsF != NULL);
   assert(varsR != NULL);
   assert(nvarsM1 >= 0 && nvarsM1 <= nvars - ntightened);
   assert(nvarsM2 >= 0 && nvarsM2 <= nvars - ntightened);
   assert(nvarsF >= 0 && nvarsF <= nvars - ntightened);
   assert(nvarsR >= 0 && nvarsR <= nvars - ntightened);
   assert(nvarsM1 + nvarsM2 + nvarsF + nvarsR == nvars  - ntightened);
   assert(alpha0 >= 0);
   assert(liftcoefs != NULL);
   assert(cutact != NULL);
   assert(liftrhs != NULL);

   /* allocates temporary memory */
   minweightssize = nvarsM1 + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &minweights, minweightssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeys, nvarsM1) );

   /* initializes data structures */
   BMSclearMemoryArray(liftcoefs, nvars);
   *cutact = 0.0;

   /* sets lifting coefficient of variables in M1, sorts variables in M1 such that a_1 <= a_2 <= ... <= a_|M1|
    * and calculates activity of the current valid inequality
    */
   for( j = 0; j < nvarsM1; j++ )
   {
      assert(liftcoefs[varsM1[j]] == 0);
      liftcoefs[varsM1[j]] = 1;
      sortkeys[j] = (SCIP_Real) (weights[varsM1[j]]);
      (*cutact) += solvals[varsM1[j]];
   }

   SCIPsortRealInt(sortkeys, varsM1, nvarsM1);

   /* initializes (i = 1) the minweight table, defined as: minweights_i[w] =
    *   min   sum_{j in M_1} a_j x_j + sum_{k=1}^{i-1} a_{j_k}     x_{j_k}
    *   s.t.  sum_{j in M_1}     x_j + sum_{k=1}^{i-1} alpha_{j_k} x_{j_k} >= w
    *                                    x_j in {0,1} for j in M_1 & {j_i,...,j_i-1},
    * for i = 1,...,t with t = |N\M1| and w = 0,...,|M1| + sum_{k=1}^{i-1} alpha_{j_k};
    */
   minweights[0] = 0;
   for( w = 1; w <= nvarsM1; w++ )
      minweights[w] = minweights[w-1] + weights[varsM1[w-1]];
   minweightslen = nvarsM1 + 1;

   /* gets sum of weights of variables fixed to one, i.e. sum of weights of variables in M_2 */
   fixedonesweight = 0;
   for( j = 0; j < nvarsM2; j++ )
      fixedonesweight += weights[varsM2[j]];
   assert(fixedonesweight >= 0);

   /* initializes right hand side of lifted valid inequality */
   *liftrhs = alpha0;

   /* sequentially up-lifts all variables in F: */
   for( j = 0; j < nvarsF; j++ )
   {
      SCIP_Longint weight;
      int liftvar;
      int liftcoef;
      int z;

      liftvar = varsF[j];
      weight = weights[liftvar];
      assert(liftvar >= 0 && liftvar < nvars);
      assert(SCIPisFeasGT(scip, solvals[liftvar], 0.0));
      assert(weight > 0);

      /* knapsack problem is infeasible:
       *   sets z = 0
       */
      if( capacity - fixedonesweight - weight < 0 )
      {
         z = 0;
      }
      /* knapsack problem is feasible:
       *   sets z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - fixedonesweight - a_{j_i}  } = liftrhs,
       *   if minweights_i[liftrhs] <= a_0 - fixedonesweight - a_{j_i}
       */
      else if( minweights[*liftrhs] <= capacity - fixedonesweight - weight )
      {
         z = *liftrhs;
      }
      /* knapsack problem is feasible:
       *   uses binary search to find z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - fixedonesweight - a_{j_i} }
       */
      else
      {
         int left;
         int right;
         int middle;

         assert((*liftrhs) + 1 >= minweightslen || minweights[(*liftrhs) + 1] > capacity - fixedonesweight - weight);
         left = 0;
         right = (*liftrhs) + 1;
         while( left < right - 1 )
         {
            middle = (left + right) / 2;
            assert(0 <= middle && middle < minweightslen);
            if( minweights[middle] <= capacity - fixedonesweight - weight )
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(0 <= left && left < minweightslen);
         assert(minweights[left] <= capacity - fixedonesweight - weight );
         assert(left == minweightslen - 1 || minweights[left+1] > capacity - fixedonesweight - weight);

         /* now z = left */
         z = left;
         assert(z <= *liftrhs);
      }

      /* calculates lifting coefficients alpha_{j_i} = liftrhs - z */
      liftcoef = (*liftrhs) - z;
      liftcoefs[liftvar] = liftcoef;
      assert(liftcoef >= 0 && liftcoef <= (*liftrhs) + 1);

      /* minweight table and activity of current valid inequality will not change, if alpha_{j_i} = 0 */
      if( liftcoef == 0 )
         continue;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];

      /* enlarges current minweight table:
       *  from minweightlen = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 entries
       *  to                  |M1| + sum_{k=1}^{i  } alpha_{j_k} + 1 entries
       * and sets minweights_i[w] = infinity for
       *  w = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 , ... , |M1| + sum_{k=1}^{i} alpha_{j_k}
       */
      SCIP_CALL( enlargeMinweights(scip, &minweights, &minweightslen, &minweightssize, minweightslen + liftcoef) );

      /* updates minweight table: minweight_i+1[w] =
       *   min{ minweights_i[w], a_{j_i}},                                 if w <  alpha_j_i
       *   min{ minweights_i[w], minweights_i[w - alpha_j_i] + a_j_i},     if w >= alpha_j_i
       */
      for( w = minweightslen - 1; w >= 0; w-- )
      {
         SCIP_Longint min;
         if( w < liftcoef )
         {
            min = MIN(minweights[w], weight);
            minweights[w] = min;
         }
         else
         {
            assert(w >= liftcoef);
            min = MIN(minweights[w], minweights[w - liftcoef] + weight);
            minweights[w] = min;
         }
      }
   }
   assert(minweights[0] == 0);

   /* sequentially down-lifts all variables in M_2: */
   for( j = 0; j < nvarsM2; j++ )
   {
      SCIP_Longint weight;
      int liftvar;
      int liftcoef;
      int left;
      int right;
      int middle;
      int z;

      liftvar = varsM2[j];
      weight = weights[liftvar];
      assert(SCIPisFeasEQ(scip, solvals[liftvar], 1.0));
      assert(liftvar >= 0 && liftvar < nvars);
      assert(weight > 0);

      /* uses binary search to find
       *   z = max { w : 0 <= w <= |M_1| + sum_{k=1}^{i-1} alpha_{j_k}, minweights_[w] <= a_0 - fixedonesweight + a_{j_i}}
       */
      left = 0;
      right = minweightslen;
      while( left < right - 1 )
      {
         middle = (left + right) / 2;
         assert(0 <= middle && middle < minweightslen);
         if( minweights[middle] <= capacity - fixedonesweight + weight )
            left = middle;
         else
            right = middle;
      }
      assert(left == right - 1);
      assert(0 <= left && left < minweightslen);
      assert(minweights[left] <= capacity - fixedonesweight + weight );
      assert(left == minweightslen - 1 || minweights[left+1] > capacity - fixedonesweight + weight);

      /* now z = left */
      z = left;
      assert(z >= *liftrhs);

      /* calculates lifting coefficients alpha_{j_i} = z - liftrhs */
      liftcoef = z - (*liftrhs);
      liftcoefs[liftvar] = liftcoef;
      assert(liftcoef >= 0);

      /* updates sum of weights of variables fixed to one */
      fixedonesweight -= weight;

      /* updates right-hand side of current valid inequality */
      (*liftrhs) += liftcoef;
      assert(*liftrhs >= alpha0);

      /* minweight table and activity of current valid inequality will not change, if alpha_{j_i} = 0 */
      if( liftcoef == 0 )
         continue;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];

      /* enlarges current minweight table:
       *  from minweightlen = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 entries
       *  to                  |M1| + sum_{k=1}^{i  } alpha_{j_k} + 1 entries
       * and sets minweights_i[w] = infinity for
       *  w = |M1| + sum_{k=1}^{i-1} alpha_{j_k} + 1 , ... , |M1| + sum_{k=1}^{i} alpha_{j_k}
       */
      SCIP_CALL( enlargeMinweights(scip, &minweights, &minweightslen, &minweightssize, minweightslen + liftcoef) );

      /* updates minweight table: minweight_i+1[w] =
       *   min{ minweights_i[w], a_{j_i}},                                 if w <  alpha_j_i
       *   min{ minweights_i[w], minweights_i[w - alpha_j_i] + a_j_i},     if w >= alpha_j_i
       */
      for( w = minweightslen - 1; w >= 0; w-- )
      {
         SCIP_Longint min;
         if( w < liftcoef )
         {
            min = MIN(minweights[w], weight);
            minweights[w] = min;
         }
         else
         {
            assert(w >= liftcoef);
            min = MIN(minweights[w], minweights[w - liftcoef] + weight);
            minweights[w] = min;
         }
      }
   }
   assert(fixedonesweight == 0);
   assert(*liftrhs >= alpha0);

   /* sequentially up-lifts all variables in R: */
   for( j = 0; j < nvarsR; j++ )
   {
      SCIP_Longint weight;
      int liftvar;
      int liftcoef;
      int z;

      liftvar = varsR[j];
      weight = weights[liftvar];
      assert(liftvar >= 0 && liftvar < nvars);
      assert(SCIPisFeasEQ(scip, solvals[liftvar], 0.0));
      assert(weight > 0);
      assert(capacity - weight >= 0);
      assert((*liftrhs) + 1 >= minweightslen || minweights[(*liftrhs) + 1] > capacity - weight);

      /* sets z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - a_{j_i} } = liftrhs,
       * if minweights_i[liftrhs] <= a_0 - a_{j_i}
       */
      if( minweights[*liftrhs] <= capacity - weight )
      {
         z = *liftrhs;
      }
      /* uses binary search to find z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - a_{j_i} }
       */
      else
      {
         int left;
         int right;
         int middle;

         left = 0;
         right = (*liftrhs) + 1;
         while( left < right - 1)
         {
            middle = (left + right) / 2;
            assert(0 <= middle && middle < minweightslen);
            if( minweights[middle] <= capacity - weight )
               left = middle;
            else
               right = middle;
         }
         assert(left == right - 1);
         assert(0 <= left && left < minweightslen);
         assert(minweights[left] <= capacity - weight );
         assert(left == minweightslen - 1 || minweights[left+1] > capacity - weight);

         /* now z = left */
         z = left;
         assert(z <= *liftrhs);
      }

      /* calculates lifting coefficients alpha_{j_i} = liftrhs - z */
      liftcoef = (*liftrhs) - z;
      liftcoefs[liftvar] = liftcoef;
      assert(liftcoef >= 0 && liftcoef <= *liftrhs);

      /* minweight table and activity of current valid inequality will not change, if alpha_{j_i} = 0 */
      if( liftcoef == 0 )
         continue;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];

      /* updates minweight table: minweight_i+1[w] =
       *   min{ minweight_i[w], a_{j_i}},                                if w <  alpha_j_i
       *   min{ minweight_i[w], minweight_i[w - alpha_j_i] + a_j_i},     if w >= alpha_j_i
       */
      for( w = *liftrhs; w >= 0; w-- )
      {
         SCIP_Longint min;
         if( w < liftcoef )
         {
            min = MIN(minweights[w], weight);
            minweights[w] = min;
         }
         else
         {
            assert(w >= liftcoef);
            min = MIN(minweights[w], minweights[w - liftcoef] + weight);
            minweights[w] = min;
         }
      }
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeys);
   SCIPfreeBufferArray(scip, &minweights);

   return SCIP_OKAY;
}

/** adds two minweight values in a safe way, i.e,, ensures no overflow */
static
SCIP_Longint safeAddMinweightsGUB(
   SCIP_Longint          val1,               /**< first value to add */
   SCIP_Longint          val2                /**< second value to add */
   )
{
   assert(val1 >= 0);
   assert(val2 >= 0);

   if( val1 >= SCIP_LONGINT_MAX || val2 >= SCIP_LONGINT_MAX )
      return SCIP_LONGINT_MAX;
   else
   {
      assert(val1 <= SCIP_LONGINT_MAX - val2);
      return (val1 + val2);
   }
}

/** computes minweights table for lifting with GUBs by combining unfished and fished tables */
static
void computeMinweightsGUB(
   SCIP_Longint*         minweights,         /**< minweight table to compute */
   SCIP_Longint*         finished,           /**< given finished table */
   SCIP_Longint*         unfinished,         /**< given unfinished table */
   int                   minweightslen       /**< length of minweight, finished, and unfinished tables */
   )
{
   int w1;
   int w2;

   /* minweights_i[w] = min{finished_i[w1] + unfinished_i[w2] : w1>=0, w2>=0, w1+w2=w};
    * note that finished and unfished arrays sorted by non-decreasing weight
    */

   /* initialize minweight with w2 = 0 */
   w2 = 0;
   assert(unfinished[w2] == 0);
   for( w1 = 0; w1 < minweightslen; w1++ )
      minweights[w1] = finished[w1];

   /* consider w2 = 1, ..., minweightslen-1 */
   for( w2 = 1; w2 < minweightslen; w2++ )
   {
      if( unfinished[w2] >= SCIP_LONGINT_MAX )
         break;

      for( w1 = 0; w1 < minweightslen - w2; w1++ )
      {
         SCIP_Longint temp;

	 temp = safeAddMinweightsGUB(finished[w1], unfinished[w2]);
	 if( temp <= minweights[w1+w2] )
	    minweights[w1+w2] = temp;
      }
   }
}

/** lifts given inequality
 *    sum_{j in C_1} x_j <= alpha_0
 *  valid for
 *    S^0 = { x in {0,1}^|C_1| : sum_{j in C_1} a_j x_j <= a_0 - sum_{j in C_2} a_j;
 *                               sum_{j in Q_i} x_j <= 1, forall i in I }
 *  to a valid inequality
 *    sum_{j in C_1} x_j + sum_{j in F} alpha_j x_j + sum_{j in C_2} alpha_j x_j + sum_{j in R} alpha_j x_j
 *    <= alpha_0 + sum_{j in C_2} alpha_j
 *  for
 *    S = { x in {0,1}^|N| : sum_{j in N} a_j x_j <= a_0; sum_{j in Q_i} x_j <= 1, forall i in I };
 *  uses sequential up-lifting   for the variables in GUB constraints in gubconsGFC1,
 *       sequential down-lifting for the variables in GUB constraints in gubconsGC2, and
 *       sequential up-lifting   for the variabels in GUB constraints in gubconsGR.
 */
static
SCIP_RETCODE sequentialUpAndDownLiftingGUB(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   ngubconscapexceed,  /**< number of GUBs with only capacity exceeding variables */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all knapsack variables */
   int*                  gubconsGC1,         /**< GUBs in GC1(GNC1+GOC1) */
   int*                  gubconsGC2,         /**< GUBs in GC2 */
   int*                  gubconsGFC1,        /**< GUBs in GFC1(GNC1+GF) */
   int*                  gubconsGR,          /**< GUBs in GR */
   int                   ngubconsGC1,        /**< number of GUBs in GC1(GNC1+GOC1) */
   int                   ngubconsGC2,        /**< number of GUBs in GC2 */
   int                   ngubconsGFC1,       /**< number of GUBs in GFC1(GNC1+GF) */
   int                   ngubconsGR,         /**< number of GUBs in GR */
   int                   alpha0,             /**< rights hand side of given valid inequality */
   int*                  liftcoefs,          /**< pointer to store lifting coefficient of vars in knapsack constraint */
   SCIP_Real*            cutact,             /**< pointer to store activity of lifted valid inequality */
   int*                  liftrhs,            /**< pointer to store right hand side of the lifted valid inequality */
   int                   maxgubvarssize      /**< maximal size of GUB constraints */
   )
{
   SCIP_Longint* minweights;
   SCIP_Longint* finished;
   SCIP_Longint* unfinished;
   int* gubconsGOC1;
   int* gubconsGNC1;
   int* liftgubvars;
   SCIP_Longint fixedonesweight;
   SCIP_Longint weight;
   SCIP_Longint weightdiff1;
   SCIP_Longint weightdiff2;
   SCIP_Longint min;
   int minweightssize;
   int minweightslen;
   int nvars;
   int varidx;
   int liftgubconsidx;
   int liftvar;
   int sumliftcoef;
   int liftcoef;
   int ngubconsGOC1;
   int ngubconsGNC1;
   int left;
   int right;
   int middle;
   int nliftgubvars;
   int tmplen;
   int tmpsize;
   int j;
   int k;
   int w;
   int z;
#ifndef NDEBUG
   int ngubconss;
   int nliftgubC1;

   assert(gubset != NULL);
   ngubconss = gubset->ngubconss;
#else
   assert(gubset != NULL);
#endif

   nvars = gubset->nvars;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(gubconsGC1 != NULL);
   assert(gubconsGC2 != NULL);
   assert(gubconsGFC1 != NULL);
   assert(gubconsGR != NULL);
   assert(ngubconsGC1 >= 0 && ngubconsGC1 <= ngubconss - ngubconscapexceed);
   assert(ngubconsGC2 >= 0 && ngubconsGC2 <= ngubconss - ngubconscapexceed);
   assert(ngubconsGFC1 >= 0 && ngubconsGFC1 <= ngubconss - ngubconscapexceed);
   assert(ngubconsGR >= 0 && ngubconsGR <= ngubconss - ngubconscapexceed);
   assert(alpha0 >= 0);
   assert(liftcoefs != NULL);
   assert(cutact != NULL);
   assert(liftrhs != NULL);

   minweightssize = ngubconsGC1+1;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &liftgubvars, maxgubvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gubconsGOC1, ngubconsGC1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gubconsGNC1, ngubconsGC1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minweights, minweightssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &finished, minweightssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &unfinished, minweightssize) );

   /* initializes data structures */
   BMSclearMemoryArray(liftcoefs, nvars);
   *cutact = 0.0;

   /* gets GOC1 and GNC1 GUBs, sets lifting coefficient of variables in C1 and calculates activity of the current
    * valid inequality
    */
   ngubconsGOC1 = 0;
   ngubconsGNC1 = 0;
   for( j = 0; j < ngubconsGC1; j++ )
   {
      if( gubset->gubconsstatus[gubconsGC1[j]] == GUBCONSSTATUS_BELONGSTOSET_GOC1 )
      {
         gubconsGOC1[ngubconsGOC1] = gubconsGC1[j];
         ngubconsGOC1++;
      }
      else
      {
         assert(gubset->gubconsstatus[gubconsGC1[j]] == GUBCONSSTATUS_BELONGSTOSET_GNC1);
         gubconsGNC1[ngubconsGNC1] = gubconsGC1[j];
         ngubconsGNC1++;
      }
      for( k = 0; k < gubset->gubconss[gubconsGC1[j]]->ngubvars
              && gubset->gubconss[gubconsGC1[j]]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_C1; k++ )
      {
         varidx = gubset->gubconss[gubconsGC1[j]]->gubvars[k];
         assert(varidx >= 0 && varidx < nvars);
         assert(liftcoefs[varidx] == 0);

         liftcoefs[varidx] = 1;
         (*cutact) += solvals[varidx];
      }
      assert(k >= 1);
   }
   assert(ngubconsGOC1 + ngubconsGFC1 + ngubconsGC2 + ngubconsGR == ngubconss - ngubconscapexceed);
   assert(ngubconsGOC1 + ngubconsGNC1 == ngubconsGC1);

   /* initialize the minweight tables, defined as: for i = 1,...,m with m = |I| and w = 0,...,|gubconsGC1|;
    * - finished_i[w] =
    *   min   sum_{k = 1,2,...,i-1} sum_{j in Q_k} a_j x_j
    *   s.t.  sum_{k = 1,2,...,i-1} sum_{j in Q_k} alpha_j x_j  >= w
    *                               sum_{j in Q_k} x_j <= 1
    *                               x_j in {0,1} forall j in Q_k forall k = 1,2,...,i-1,
    * - unfinished_i[w] =
    *   min   sum_{k = i+1,...,m} sum_{j in Q_k && j in C1} a_j x_j
    *   s.t.  sum_{k = i+1,...,m} sum_{j in Q_k && j in C1} x_j  >= w
    *                             sum_{j in Q_k} x_j <= 1
    *                             x_j in {0,1} forall j in Q_k forall k = 1,2,...,i-1,
    * - minweights_i[w] = min{finished_i[w1] + unfinished_i[w2] : w1>=0, w2>=0, w1+w2=w};
    */

   /* initialize finished table; note that variables in GOC1 GUBs (includes C1 and capacity exceeding variables)
    * are sorted s.t. C1 variables come first and are sorted by non-decreasing weight.
    * GUBs in the group GCI are sorted by non-decreasing min{ a_k : k in GC1_j } where min{ a_k : k in GC1_j } always
    * comes from the first variable in the GUB
    */
   assert(ngubconsGOC1 <= ngubconsGC1);
   finished[0] = 0;
   for( w = 1; w <= ngubconsGOC1; w++ )
   {
      liftgubconsidx = gubconsGOC1[w-1];

      assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GOC1);
      assert(gubset->gubconss[liftgubconsidx]->gubvarsstatus[0] == GUBVARSTATUS_BELONGSTOSET_C1);

      varidx = gubset->gubconss[liftgubconsidx]->gubvars[0];

      assert(varidx >= 0 && varidx < nvars);
      assert(liftcoefs[varidx] == 1);

      min = weights[varidx];
      finished[w] = finished[w-1] + min;

#ifndef NDEBUG
      for( k = 1; k < gubset->gubconss[liftgubconsidx]->ngubvars
              && gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_C1; k++ )
      {
         varidx = gubset->gubconss[liftgubconsidx]->gubvars[k];
         assert(varidx >= 0 && varidx < nvars);
         assert(liftcoefs[varidx] == 1);
         assert(weights[varidx] >= min);
      }
#endif
   }
   for( w = ngubconsGOC1+1; w <= ngubconsGC1; w++ )
      finished[w] = SCIP_LONGINT_MAX;

   /* initialize unfinished table; note that variables in GNC1 GUBs
    * are sorted s.t. C1 variables come first and are sorted by non-decreasing weight.
    * GUBs in the group GCI are sorted by non-decreasing min{ a_k : k in GC1_j } where min{ a_k : k in GC1_j } always
    * comes from the first variable in the GUB
    */
   assert(ngubconsGNC1 <= ngubconsGC1);
   unfinished[0] = 0;
   for( w = 1; w <= ngubconsGNC1; w++ )
   {
      liftgubconsidx = gubconsGNC1[w-1];

      assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1);
      assert(gubset->gubconss[liftgubconsidx]->gubvarsstatus[0] == GUBVARSTATUS_BELONGSTOSET_C1);

      varidx = gubset->gubconss[liftgubconsidx]->gubvars[0];

      assert(varidx >= 0 && varidx < nvars);
      assert(liftcoefs[varidx] == 1);

      min = weights[varidx];
      unfinished[w] = unfinished[w-1] + min;

#ifndef NDEBUG
      for( k = 1; k < gubset->gubconss[liftgubconsidx]->ngubvars
              && gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_C1; k++ )
      {
         varidx = gubset->gubconss[liftgubconsidx]->gubvars[k];
         assert(varidx >= 0 && varidx < nvars);
         assert(liftcoefs[varidx] == 1);
         assert(weights[varidx] >= min );
      }
#endif
   }
   for( w = ngubconsGNC1 + 1; w <= ngubconsGC1; w++ )
      unfinished[w] = SCIP_LONGINT_MAX;

   /* initialize minweights table; note that variables in GC1 GUBs
    * are sorted s.t. C1 variables come first and are sorted by non-decreasing weight.
    * we can directly initialize minweights instead of computing it from finished and unfinished (which would be more time
    * consuming) because is it has to be build using weights from C1 only.
    */
   assert(ngubconsGOC1 + ngubconsGNC1 == ngubconsGC1);
   minweights[0] = 0;
   for( w = 1; w <= ngubconsGC1; w++ )
   {
      liftgubconsidx = gubconsGC1[w-1];

      assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GOC1
          || gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1);
      assert(gubset->gubconss[liftgubconsidx]->gubvarsstatus[0] == GUBVARSTATUS_BELONGSTOSET_C1);

      varidx = gubset->gubconss[liftgubconsidx]->gubvars[0];

      assert(varidx >= 0 && varidx < nvars);
      assert(liftcoefs[varidx] == 1);

      min = weights[varidx];
      minweights[w] = minweights[w-1] + min;

#ifndef NDEBUG
      for( k = 1; k < gubset->gubconss[liftgubconsidx]->ngubvars
              && gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_C1; k++ )
      {
         varidx = gubset->gubconss[liftgubconsidx]->gubvars[k];
         assert(varidx >= 0 && varidx < nvars);
         assert(liftcoefs[varidx] == 1);
         assert(weights[varidx] >= min);
      }
#endif
   }
   minweightslen = ngubconsGC1 + 1;

   /* gets sum of weights of variables fixed to one, i.e. sum of weights of C2 variables GC2 GUBs */
   fixedonesweight = 0;
   for( j = 0; j < ngubconsGC2; j++ )
   {
      varidx = gubset->gubconss[gubconsGC2[j]]->gubvars[0];

      assert(gubset->gubconss[gubconsGC2[j]]->ngubvars == 1);
      assert(varidx >= 0 && varidx < nvars);
      assert(gubset->gubconss[gubconsGC2[j]]->gubvarsstatus[0] == GUBVARSTATUS_BELONGSTOSET_C2);

      fixedonesweight += weights[varidx];
   }
   assert(fixedonesweight >= 0);

   /* initializes right hand side of lifted valid inequality */
   *liftrhs = alpha0;

   /* sequentially up-lifts all variables in GFC1 GUBs */
   for( j = 0; j < ngubconsGFC1; j++ )
   {
      liftgubconsidx = gubconsGFC1[j];
      assert(liftgubconsidx >= 0 && liftgubconsidx < ngubconss);

      /* GNC1 GUB: update unfinished table (remove current GUB, i.e., remove min weight of C1 vars in GUB) and
       * compute minweight table via updated unfinished table and aleady upto date finished table;
       */
      k = 0;
      if( gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1 )
      {
	 assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1);
         assert(gubset->gubconss[liftgubconsidx]->gubvarsstatus[0] == GUBVARSTATUS_BELONGSTOSET_C1);
         assert(ngubconsGNC1 > 0);

         /* get number of C1 variables of current GNC1 GUB and put them into array of variables in GUB that
          * are considered for the lifting, i.e., not capacity exceeding
          */
         for( ; k < gubset->gubconss[liftgubconsidx]->ngubvars
		&& gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_C1; k++ )
            liftgubvars[k] = gubset->gubconss[liftgubconsidx]->gubvars[k];
         assert(k >= 1);

         /* update unfinished table by removing current GNC1 GUB, i.e, remove C1 variable with minimal weight
	  * unfinished[w] = MAX{unfinished[w], unfinished[w+1] - weight}, "weight" is the minimal weight of current GUB
	  */
         weight = weights[liftgubvars[0]];

	 weightdiff2 = unfinished[ngubconsGNC1] - weight;
	 unfinished[ngubconsGNC1] = SCIP_LONGINT_MAX;
         for( w = ngubconsGNC1-1; w >= 1; w-- )
         {
	    weightdiff1 = weightdiff2;
	    weightdiff2 = unfinished[w] - weight;

            if( unfinished[w] < weightdiff1 )
	       unfinished[w] = weightdiff1;
	    else
	       break;
         }
         ngubconsGNC1--;

         /* computes minweights table by combining unfished and fished tables */
         computeMinweightsGUB(minweights, finished, unfinished, minweightslen);
         assert(minweights[0] == 0);
      }
      /* GF GUB: no update of unfinished table (and minweight table) required because GF GUBs have no C1 variables and
       * are therefore not in the unfinished table
       */
      else
	 assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GF);

#ifndef NDEBUG
      nliftgubC1 = k;
#endif
      nliftgubvars = k;
      sumliftcoef = 0;

      /* compute lifting coefficient of F and R variables in GNC1 and GF GUBs (C1 vars have already liftcoef 1) */
      for( ; k < gubset->gubconss[liftgubconsidx]->ngubvars; k++ )
      {
         if( gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_F
             || gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_R )
         {
            liftvar = gubset->gubconss[liftgubconsidx]->gubvars[k];
            weight = weights[liftvar];
            assert(weight > 0);
            assert(liftvar >= 0 && liftvar < nvars);
	    assert(capacity - weight >= 0);

            /* put variable into array of variables in GUB that are considered for the lifting,
             *  i.e., not capacity exceeding
             */
            liftgubvars[nliftgubvars] = liftvar;
            nliftgubvars++;

	    /* knapsack problem is infeasible:
             * sets z = 0
             */
            if( capacity - fixedonesweight - weight < 0 )
            {
               z = 0;
            }
            /* knapsack problem is feasible:
             *   sets z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - fixedonesweight - a_{j_i}  } = liftrhs,
             *   if minweights_i[liftrhs] <= a_0 - fixedonesweight - a_{j_i}
             */
            else if( minweights[*liftrhs] <= capacity - fixedonesweight - weight )
            {
               z = *liftrhs;
            }
            /* knapsack problem is feasible:
             *   binary search to find z = max {w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - fixedonesweight - a_{j_i}}
             */
            else
            {
               assert((*liftrhs) + 1 >= minweightslen || minweights[(*liftrhs) + 1] > capacity - fixedonesweight - weight);
               left = 0;
               right = (*liftrhs) + 1;
               while( left < right - 1 )
               {
                  middle = (left + right) / 2;
                  assert(0 <= middle && middle < minweightslen);
                  if( minweights[middle] <= capacity - fixedonesweight - weight )
                     left = middle;
                  else
                     right = middle;
               }
               assert(left == right - 1);
               assert(0 <= left && left < minweightslen);
               assert(minweights[left] <= capacity - fixedonesweight - weight);
               assert(left == minweightslen - 1 || minweights[left+1] > capacity - fixedonesweight - weight);

               /* now z = left */
               z = left;
               assert(z <= *liftrhs);
            }

            /* calculates lifting coefficients alpha_{j_i} = liftrhs - z */
            liftcoef = (*liftrhs) - z;
            liftcoefs[liftvar] = liftcoef;
            assert(liftcoef >= 0 && liftcoef <= (*liftrhs) + 1);

	    /* updates activity of current valid inequality */
            (*cutact) += liftcoef * solvals[liftvar];

            /* updates sum of all lifting coefficients in GUB */
            sumliftcoef += liftcoefs[liftvar];
	 }
         else
            assert(gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_CAPACITYEXCEEDED);
      }
      /* at least one variable is in F or R (j = number of C1 variables in current GUB) */
      assert(nliftgubvars > nliftgubC1);

      /* activity of current valid inequality will not change if (sum of alpha_{j_i} in GUB) = 0
       * and finished and minweight table can be updated easily as only C1 variables need to be considered;
       * not needed for GF GUBs
       */
      if( sumliftcoef == 0 )
      {
	 if( gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GNC1 )
	 {
            weight = weights[liftgubvars[0]];
            /* update finished table and minweights table by applying special case of
             * finished[w] = MIN{finished[w], finished[w-1] + weight}, "weight" is the minimal weight of current GUB
	     * minweights[w] = MIN{minweights[w], minweights[w-1] + weight}, "weight" is the minimal weight of current GUB
             */
            for( w = minweightslen-1; w >= 1; w-- )
            {
               SCIP_Longint tmpval;

               tmpval = safeAddMinweightsGUB(finished[w-1], weight);
               finished[w] = MIN(finished[w], tmpval);

               tmpval = safeAddMinweightsGUB(minweights[w-1], weight);
               minweights[w] = MIN(minweights[w], tmpval);
            }
         }
         else
	    assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GF);

         continue;
      }

      /* enlarges current minweights tables(finished, unfinished, minweights):
       *  from minweightlen = |gubconsGC1| + sum_{k=1,2,...,i-1}sum_{j in Q_k} alpha_j + 1 entries
       *  to                  |gubconsGC1| + sum_{k=1,2,...,i  }sum_{j in Q_k} alpha_j + 1 entries
       *  and sets minweights_i[w] = infinity for
       *  w = |gubconsGC1| + sum_{k=1,2,..,i-1}sum_{j in Q_k} alpha_j+1,..,|C1| + sum_{k=1,2,..,i}sum_{j in Q_k} alpha_j
       */
      tmplen = minweightslen; /* will be updated in enlargeMinweights() */
      tmpsize = minweightssize;
      SCIP_CALL( enlargeMinweights(scip, &unfinished, &tmplen, &tmpsize, tmplen + sumliftcoef) );
      tmplen = minweightslen;
      tmpsize = minweightssize;
      SCIP_CALL( enlargeMinweights(scip, &finished, &tmplen, &tmpsize, tmplen + sumliftcoef) );
      SCIP_CALL( enlargeMinweights(scip, &minweights, &minweightslen, &minweightssize, minweightslen + sumliftcoef) );

      /* update finished table and minweight table;
       * note that instead of computing minweight table from updated finished and updated unfinished table again
       * (for the lifting coefficient, we had to update unfinished table and compute minweight table), we here
       * only need to update the minweight table and the updated finished in the same way (i.e., computing for minweight
       * not needed because only finished table changed at this point and the change was "adding" one weight)
       *
       * update formular for minweight table is: minweight_i+1[w] =
       *   min{ minweights_i[w], min{ minweights_i[w - alpha_k]^{+} + a_k : k in GUB_j_i } }
       * formular for finished table has the same pattern.
       */
      for( w = minweightslen-1; w >= 0; w-- )
      {
         SCIP_Longint minminweight;
         SCIP_Longint minfinished;

         for( k = 0; k < nliftgubvars; k++ )
	 {
	    liftcoef = liftcoefs[liftgubvars[k]];
	    weight = weights[liftgubvars[k]];

            if( w < liftcoef )
            {
	       minfinished = MIN(finished[w], weight);
	       minminweight = MIN(minweights[w], weight);

               finished[w] = minfinished;
               minweights[w] = minminweight;
            }
            else
            {
               SCIP_Longint tmpval;

               assert(w >= liftcoef);

               tmpval = safeAddMinweightsGUB(finished[w-liftcoef], weight);
               minfinished = MIN(finished[w], tmpval);

               tmpval = safeAddMinweightsGUB(minweights[w-liftcoef], weight);
               minminweight = MIN(minweights[w], tmpval);

               finished[w] = minfinished;
               minweights[w] = minminweight;
            }
	 }
      }
      assert(minweights[0] == 0);
   }
   assert(ngubconsGNC1 == 0);

   /* note: now the unfinished table no longer exists, i.e., it is "0, MAX, MAX, ..." and minweight equals to finished;
    * therefore, only work with minweight table from here on
    */

   /* sequentially down-lifts C2 variables contained in trivial GC2 GUBs */
   for( j = 0; j < ngubconsGC2; j++ )
   {
      liftgubconsidx = gubconsGC2[j];

      assert(liftgubconsidx >=0 && liftgubconsidx < ngubconss);
      assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GC2);
      assert(gubset->gubconss[liftgubconsidx]->ngubvars == 1);
      assert(gubset->gubconss[liftgubconsidx]->gubvarsstatus[0] == GUBVARSTATUS_BELONGSTOSET_C2);

      liftvar = gubset->gubconss[liftgubconsidx]->gubvars[0]; /* C2 GUBs contain only one variable */
      weight = weights[liftvar];

      assert(liftvar >= 0 && liftvar < nvars);
      assert(SCIPisFeasEQ(scip, solvals[liftvar], 1.0));
      assert(weight > 0);

      /* uses binary search to find
       *   z = max { w : 0 <= w <= |C_1| + sum_{k=1}^{i-1} alpha_{j_k}, minweights_[w] <= a_0 - fixedonesweight + a_{j_i}}
       */
      left = 0;
      right = minweightslen;
      while( left < right - 1 )
      {
         middle = (left + right) / 2;
         assert(0 <= middle && middle < minweightslen);
         if( minweights[middle] <= capacity - fixedonesweight + weight )
            left = middle;
         else
            right = middle;
      }
      assert(left == right - 1);
      assert(0 <= left && left < minweightslen);
      assert(minweights[left] <= capacity - fixedonesweight + weight);
      assert(left == minweightslen - 1 || minweights[left + 1] > capacity - fixedonesweight + weight);

      /* now z = left */
      z = left;
      assert(z >= *liftrhs);

      /* calculates lifting coefficients alpha_{j_i} = z - liftrhs */
      liftcoef = z - (*liftrhs);
      liftcoefs[liftvar] = liftcoef;
      assert(liftcoef >= 0);

      /* updates sum of weights of variables fixed to one */
      fixedonesweight -= weight;

      /* updates right-hand side of current valid inequality */
      (*liftrhs) += liftcoef;
      assert(*liftrhs >= alpha0);

      /* minweight table and activity of current valid inequality will not change, if alpha_{j_i} = 0 */
      if( liftcoef == 0 )
	 continue;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];

      /* enlarges current minweight table:
       *  from minweightlen = |gubconsGC1| + sum_{k=1,2,...,i-1}sum_{j in Q_k} alpha_j + 1 entries
       *  to                  |gubconsGC1| + sum_{k=1,2,...,i  }sum_{j in Q_k} alpha_j + 1 entries
       * and sets minweights_i[w] = infinity for
       *  w = |C1| + sum_{k=1,2,...,i-1}sum_{j in Q_k} alpha_j + 1 , ... , |C1| + sum_{k=1,2,...,i}sum_{j in Q_k} alpha_j
       */
      SCIP_CALL( enlargeMinweights(scip, &minweights, &minweightslen, &minweightssize, minweightslen + liftcoef) );

      /* updates minweight table: minweight_i+1[w] =
       *  min{ minweights_i[w], a_{j_i}},                                 if w <  alpha_j_i
       *  min{ minweights_i[w], minweights_i[w - alpha_j_i] + a_j_i},     if w >= alpha_j_i
       */
      for( w = minweightslen - 1; w >= 0; w-- )
      {
         if( w < liftcoef )
         {
            min = MIN(minweights[w], weight);
            minweights[w] = min;
         }
         else
         {
            SCIP_Longint tmpval;

            assert(w >= liftcoef);

            tmpval = safeAddMinweightsGUB(minweights[w-liftcoef], weight);
            min = MIN(minweights[w], tmpval);
            minweights[w] = min;
         }
      }
   }
   assert(fixedonesweight == 0);
   assert(*liftrhs >= alpha0);

   /* sequentially up-lifts variables in GUB constraints in GR GUBs */
   for( j = 0; j < ngubconsGR; j++ )
   {
      liftgubconsidx = gubconsGR[j];

      assert(liftgubconsidx >=0 && liftgubconsidx < ngubconss);
      assert(gubset->gubconsstatus[liftgubconsidx] == GUBCONSSTATUS_BELONGSTOSET_GR);

      sumliftcoef = 0;
      nliftgubvars = 0;
      for( k = 0; k < gubset->gubconss[liftgubconsidx]->ngubvars; k++ )
      {
         if(gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_BELONGSTOSET_R )
         {
            liftvar = gubset->gubconss[liftgubconsidx]->gubvars[k];
            weight = weights[liftvar];
            assert(weight > 0);
            assert(liftvar >= 0 && liftvar < nvars);
	    assert(capacity - weight >= 0);
            assert((*liftrhs) + 1 >= minweightslen || minweights[(*liftrhs) + 1] > capacity - weight);

            /* put variable into array of variables in GUB that are considered for the lifting,
             *  i.e., not capacity exceeding
             */
	    liftgubvars[nliftgubvars] = liftvar;
            nliftgubvars++;

            /* sets z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - a_{j_i} } = liftrhs,
             * if minweights_i[liftrhs] <= a_0 - a_{j_i}
             */
            if( minweights[*liftrhs] <= capacity - weight )
            {
               z = *liftrhs;
            }
            /* uses binary search to find z = max { w : 0 <= w <= liftrhs, minweights_i[w] <= a_0 - a_{j_i} }
             */
            else
            {
               left = 0;
               right = (*liftrhs) + 1;
               while( left < right - 1 )
               {
                  middle = (left + right) / 2;
                  assert(0 <= middle && middle < minweightslen);
                  if( minweights[middle] <= capacity - weight )
                     left = middle;
                  else
                     right = middle;
               }
               assert(left == right - 1);
               assert(0 <= left && left < minweightslen);
               assert(minweights[left] <= capacity - weight);
               assert(left == minweightslen - 1 || minweights[left + 1] > capacity - weight);

               /* now z = left */
               z = left;
               assert(z <= *liftrhs);
            }
            /* calculates lifting coefficients alpha_{j_i} = liftrhs - z */
            liftcoef = (*liftrhs) - z;
            liftcoefs[liftvar] = liftcoef;
            assert(liftcoef >= 0 && liftcoef <= (*liftrhs) + 1);

            /* updates activity of current valid inequality */
            (*cutact) += liftcoef * solvals[liftvar];

            /* updates sum of all lifting coefficients in GUB */
            sumliftcoef += liftcoefs[liftvar];
         }
         else
            assert(gubset->gubconss[liftgubconsidx]->gubvarsstatus[k] == GUBVARSTATUS_CAPACITYEXCEEDED);
      }
      assert(nliftgubvars >= 1); /* at least one variable is in R */

      /* minweight table and activity of current valid inequality will not change if (sum of alpha_{j_i} in GUB) = 0 */
      if( sumliftcoef == 0 )
	 continue;

      /* updates minweight table: minweight_i+1[w] =
       *   min{ minweights_i[w], min{ minweights_i[w - alpha_k]^{+} + a_k : k in GUB_j_i } }
       */
      for( w = *liftrhs; w >= 0; w-- )
      {
         for( k = 0; k < nliftgubvars; k++ )
         {
            liftcoef = liftcoefs[liftgubvars[k]];
            weight = weights[liftgubvars[k]];

            if( w < liftcoef )
            {
               min = MIN(minweights[w], weight);
               minweights[w] = min;
            }
            else
            {
               SCIP_Longint tmpval;

               assert(w >= liftcoef);

               tmpval = safeAddMinweightsGUB(minweights[w-liftcoef], weight);
               min = MIN(minweights[w], tmpval);
               minweights[w] = min;
            }
         }
      }
      assert(minweights[0] == 0);
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &minweights);
   SCIPfreeBufferArray(scip, &finished);
   SCIPfreeBufferArray(scip, &unfinished);
   SCIPfreeBufferArray(scip, &liftgubvars);
   SCIPfreeBufferArray(scip, &gubconsGOC1 );
   SCIPfreeBufferArray(scip, &gubconsGNC1);

   return SCIP_OKAY;
}

/** lifts given minimal cover inequality 
 *  \f[
 *    \sum_{j \in C} x_j \leq |C| - 1 
 *  \f]
 *  valid for 
 *  \f[
 *    S^0 = \{ x \in {0,1}^{|C|} : \sum_{j \in C} a_j x_j \leq a_0 \}
 *  \f]
 *  to a valid inequality 
 *  \f[
 *    \sum_{j \in C} x_j + \sum_{j \in N \setminus C} \alpha_j x_j \leq |C| - 1
 *  \f]
 *  for
 *  \f[ 
 *    S = \{ x \in {0,1}^{|N|} : \sum_{j \in N} a_j x_j \leq a_0 \}; 
 *  \f]
 *  uses superadditive up-lifting for the variables in \f$N \setminus C\f$.
 */ 
static
SCIP_RETCODE superadditiveUpLifting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables */
   int*                  noncovervars,       /**< noncover variables */
   int                   ncovervars,         /**< number of cover variables */
   int                   nnoncovervars,      /**< number of noncover variables */
   SCIP_Longint          coverweight,        /**< weight of cover */
   SCIP_Real*            liftcoefs,          /**< pointer to store lifting coefficient of vars in knapsack constraint */
   SCIP_Real*            cutact              /**< pointer to store activity of lifted valid inequality */
   )
{
   SCIP_Longint* maxweightsums;
   SCIP_Longint* intervalends;
   SCIP_Longint* rhos;
   SCIP_Real* sortkeys;
   SCIP_Longint lambda;
   int j;
   int h;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(weights != NULL);
   assert(capacity >= 0);
   assert(solvals != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars > 0 && ncovervars <= nvars);
   assert(nnoncovervars >= 0 && nnoncovervars <= nvars - ntightened);
   assert(ncovervars + nnoncovervars == nvars - ntightened);
   assert(liftcoefs != NULL);
   assert(cutact != NULL);

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeys, ncovervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxweightsums, ncovervars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &intervalends, ncovervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhos, ncovervars) );

   /* initializes data structures */
   BMSclearMemoryArray(liftcoefs, nvars);
   *cutact = 0.0;

   /* sets lifting coefficient of variables in C, sorts variables in C such that a_1 >= a_2 >= ... >= a_|C| 
    * and calculates activity of current valid inequality 
    */
   for( j = 0; j < ncovervars; j++ )
   {
      assert(liftcoefs[covervars[j]] == 0.0);
      liftcoefs[covervars[j]] = 1.0;
      sortkeys[j] = (SCIP_Real) weights[covervars[j]];
      (*cutact) += solvals[covervars[j]];
   }
   SCIPsortDownRealInt(sortkeys, covervars, ncovervars);

   /* calculates weight excess of cover C */
   lambda = coverweight - capacity;
   assert(lambda > 0);

   /* calculates A_h for h = 0,...,|C|, I_h for h = 1,...,|C| and rho_h for h = 1,...,|C| */
   maxweightsums[0] = 0;
   for( h = 1; h <= ncovervars; h++ )
   {
      maxweightsums[h] = maxweightsums[h-1] + weights[covervars[h-1]];
      intervalends[h-1] = maxweightsums[h] - lambda;
      rhos[h-1] = MAX(0, weights[covervars[h-1]] - weights[covervars[0]] + lambda);
   }

   /* sorts variables in N\C such that a_{j_1} <= a_{j_2} <= ... <= a_{j_t} */
   for( j = 0; j < nnoncovervars; j++ )
      sortkeys[j] = (SCIP_Real) (weights[noncovervars[j]]);
   SCIPsortRealInt(sortkeys, noncovervars, nnoncovervars);

   /* calculates lifting coefficient for all variables in N\C */
   h = 0;
   for( j = 0; j < nnoncovervars; j++ )
   {
      int liftvar;
      SCIP_Longint weight;
      SCIP_Real liftcoef;

      liftvar = noncovervars[j];
      weight = weights[liftvar];

      while( intervalends[h] < weight )
         h++;

      if( h == 0 )
         liftcoef = h;
      else 
      {
         if( weight <= intervalends[h-1] + rhos[h] )
         {
            SCIP_Real tmp1;
            SCIP_Real tmp2;
            tmp1 =  (SCIP_Real) (intervalends[h-1] + rhos[h] - weight);
            tmp2 =  (SCIP_Real) rhos[1];
            liftcoef = h - ( tmp1 / tmp2 );
         }
         else
            liftcoef = h;
      }      

      /* sets lifting coefficient */
      assert(liftcoefs[liftvar] == 0.0);
      liftcoefs[liftvar] = liftcoef;

      /* updates activity of current valid inequality */
      (*cutact) += liftcoef * solvals[liftvar];
   } 

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &rhos);
   SCIPfreeBufferArray(scip, &intervalends);
   SCIPfreeBufferArray(scip, &maxweightsums);
   SCIPfreeBufferArray(scip, &sortkeys);

   return SCIP_OKAY;
}


/** separates lifted minimal cover inequalities using sequential up- and down-lifting and GUB information, if wanted, for
 *  given knapsack problem
*/
static
SCIP_RETCODE separateSequLiftedMinimalCoverInequality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< originating constraint of the knapsack problem, or NULL */
   SCIP_SEPA*            sepa,               /**< originating separator of the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  mincovervars,       /**< mincover variables */
   int*                  nonmincovervars,    /**< nonmincover variables */
   int                   nmincovervars,      /**< number of mincover variables */
   int                   nnonmincovervars,   /**< number of nonmincover variables */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   SCIP_GUBSET*          gubset,             /**< GUB set data structure, NULL if no GUB information should be used */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   int* varsC1;
   int* varsC2;
   int* varsF;
   int* varsR;
   int nvarsC1;
   int nvarsC2;
   int nvarsF;
   int nvarsR;
   SCIP_Real cutact;
   int* liftcoefs;
   int liftrhs;

   assert( cutoff != NULL );
   *cutoff = FALSE;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &varsC1, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsC2, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsF, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsR, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcoefs, nvars) );

   /* gets partition (C_1,C_2) of C, i.e. C_1 & C_2 = C and C_1 cap C_2 = emptyset, with C_1 not empty; chooses partition
    * as follows
    *   C_2 = { j in C : x*_j = 1 } and
    *   C_1 = C\C_2
    */
   getPartitionCovervars(scip, solvals, mincovervars, nmincovervars, varsC1, varsC2, &nvarsC1, &nvarsC2);
   assert(nvarsC1 + nvarsC2 == nmincovervars);
   assert(nmincovervars > 0);
   assert(nvarsC1 >= 0); /* nvarsC1 > 0 does not always hold, because relaxed knapsack conss may already be violated */

   /* changes partition (C_1,C_2) of minimal cover C, if |C1| = 1, by moving one variable from C2 to C1 */
   if( nvarsC1 < 2 && nvarsC2 > 0)
   {
      SCIP_CALL( changePartitionCovervars(scip, weights, varsC1, varsC2, &nvarsC1, &nvarsC2) );
      assert(nvarsC1 >= 1);
   }
   assert(nvarsC2 == 0 || nvarsC1 >= 1);

   /* gets partition (F,R) of N\C, i.e. F & R = N\C and F cap R = emptyset; chooses partition as follows
    *   R = { j in N\C : x*_j = 0 } and
    *   F = (N\C)\F
    */
   getPartitionNoncovervars(scip, solvals, nonmincovervars, nnonmincovervars, varsF, varsR, &nvarsF, &nvarsR);
   assert(nvarsF + nvarsR == nnonmincovervars);
   assert(nvarsC1 + nvarsC2 + nvarsF + nvarsR == nvars - ntightened);

   /* lift cuts without GUB information */
   if( gubset == NULL )
   {
      /* sorts variables in F, C_2, R according to the second level lifting sequence that will be used in the sequential
       * lifting procedure
       */
      SCIP_CALL( getLiftingSequence(scip, solvals, weights, varsF, varsC2, varsR, nvarsF, nvarsC2, nvarsR) );

      /* lifts minimal cover inequality sum_{j in C_1} x_j <= |C_1| - 1 valid for
       *
       *    S^0 = { x in {0,1}^|C_1| : sum_{j in C_1} a_j x_j <= a_0 - sum_{j in C_2} a_j }
       *
       * to a valid inequality sum_{j in C_1} x_j + sum_{j in N\C_1} alpha_j x_j <= |C_1| - 1 + sum_{j in C_2} alpha_j for
       *
       *      S = { x in {0,1}^|N|   : sum_{j in N}   a_j x_j <= a_0 },
       *
       * uses sequential up-lifting for the variables in F, sequential down-lifting for the variable in C_2 and sequential
       * up-lifting for the variables in R according to the second level lifting sequence
       */
      SCIP_CALL( sequentialUpAndDownLifting(scip, vars, nvars, ntightened, weights, capacity, solvals, varsC1, varsC2,
            varsF, varsR, nvarsC1, nvarsC2, nvarsF, nvarsR, nvarsC1 - 1, liftcoefs, &cutact, &liftrhs) );
   }
   /* lift cuts with GUB information */
   else
   {
      int* gubconsGC1;
      int* gubconsGC2;
      int* gubconsGFC1;
      int* gubconsGR;
      int ngubconsGC1;
      int ngubconsGC2;
      int ngubconsGFC1;
      int ngubconsGR;
      int ngubconss;
      int nconstightened;
      int maxgubvarssize;

      assert(nvars == gubset->nvars);

      ngubconsGC1 = 0;
      ngubconsGC2 = 0;
      ngubconsGFC1 = 0;
      ngubconsGR = 0;
      ngubconss = gubset->ngubconss;
      nconstightened = 0;
      maxgubvarssize = 0;

      /* allocates temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &gubconsGC1, ngubconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &gubconsGC2, ngubconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &gubconsGFC1, ngubconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &gubconsGR, ngubconss) );

      /* categorizies GUBs of knapsack GUB partion into GOC1, GNC1, GF, GC2, and GR and computes a lifting sequence of
       * the GUBs for the sequential GUB wise lifting procedure
       */
      SCIP_CALL( getLiftingSequenceGUB(scip, gubset, solvals, weights, varsC1, varsC2, varsF, varsR, nvarsC1,
            nvarsC2, nvarsF, nvarsR, gubconsGC1, gubconsGC2, gubconsGFC1, gubconsGR, &ngubconsGC1, &ngubconsGC2,
            &ngubconsGFC1, &ngubconsGR, &nconstightened, &maxgubvarssize) );

      /* lifts minimal cover inequality sum_{j in C_1} x_j <= |C_1| - 1 valid for
       *
       *    S^0 = { x in {0,1}^|C_1| : sum_{j in C_1} a_j x_j <= a_0 - sum_{j in C_2} a_j,
       *                               sum_{j in Q_i} x_j <= 1, forall i in I }
       *
       * to a valid inequality sum_{j in C_1} x_j + sum_{j in N\C_1} alpha_j x_j <= |C_1| - 1 + sum_{j in C_2} alpha_j for
       *
       *      S = { x in {0,1}^|N|   : sum_{j in N}   a_j x_j <= a_0, sum_{j in Q_i} x_j <= 1, forall i in I },
       *
       *  uses sequential up-lifting   for the variables in GUB constraints in gubconsGFC1,
       *       sequential down-lifting for the variables in GUB constraints in gubconsGC2, and
       *       sequential up-lifting   for the variabels in GUB constraints in gubconsGR.
       */
      SCIP_CALL( sequentialUpAndDownLiftingGUB(scip, gubset, vars, nconstightened, weights, capacity, solvals, gubconsGC1,
            gubconsGC2, gubconsGFC1, gubconsGR, ngubconsGC1, ngubconsGC2, ngubconsGFC1, ngubconsGR,
            MIN(nvarsC1 - 1, ngubconsGC1), liftcoefs, &cutact, &liftrhs, maxgubvarssize) );

      /* frees temporary memory */
      SCIPfreeBufferArray(scip, &gubconsGR);
      SCIPfreeBufferArray(scip, &gubconsGFC1);
      SCIPfreeBufferArray(scip, &gubconsGC2);
      SCIPfreeBufferArray(scip, &gubconsGC1);
   }

   /* checks, if lifting yielded a violated cut */
   if( SCIPisEfficacious(scip, (cutact - liftrhs)/sqrt((SCIP_Real)MAX(liftrhs, 1))) )
   {
      SCIP_ROW* row;
      char name[SCIP_MAXSTRLEN];
      int j;

      /* creates LP row */
      assert( cons == NULL || sepa == NULL );
      if ( cons != NULL )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_mcseq%" SCIP_LONGINT_FORMAT "", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real)liftrhs,
            cons != NULL ? SCIPconsIsLocal(cons) : FALSE, FALSE,
            cons != NULL ? SCIPconsIsRemovable(cons) : TRUE) );
      }
      else if ( sepa != NULL )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_mcseq_%" SCIP_LONGINT_FORMAT "", SCIPsepaGetName(sepa), SCIPsepaGetNCutsFound(sepa));
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, FALSE, FALSE, TRUE) );
      }
      else
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nn_mcseq_%" SCIP_LONGINT_FORMAT "", *ncuts);
         SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, FALSE, FALSE, TRUE) );
      }

      /* adds all variables in the knapsack constraint with calculated lifting coefficient to the cut */
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      assert(nvarsC1 + nvarsC2 + nvarsF + nvarsR == nvars - ntightened);
      for( j = 0; j < nvarsC1; j++ )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsC1[j]], 1.0) );
      }
      for( j = 0; j < nvarsC2; j++ )
      {
         if( liftcoefs[varsC2[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsC2[j]], (SCIP_Real)liftcoefs[varsC2[j]]) );
         }
      }
      for( j = 0; j < nvarsF; j++ )
      {
         if( liftcoefs[varsF[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsF[j]], (SCIP_Real)liftcoefs[varsF[j]]) );
         }
      }
      for( j = 0; j < nvarsR; j++ )
      {
         if( liftcoefs[varsR[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsR[j]], (SCIP_Real)liftcoefs[varsR[j]]) );
         }
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      /* checks, if cut is violated enough */
      if( SCIPisCutEfficacious(scip, sol, row) )
      {
         if( cons != NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
         (*ncuts)++;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &liftcoefs);
   SCIPfreeBufferArray(scip, &varsR);
   SCIPfreeBufferArray(scip, &varsF);
   SCIPfreeBufferArray(scip, &varsC2);
   SCIPfreeBufferArray(scip, &varsC1);

   return SCIP_OKAY;
}

/** separates lifted extended weight inequalities using sequential up- and down-lifting for given knapsack problem */
static
SCIP_RETCODE separateSequLiftedExtendedWeightInequality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_SEPA*            sepa,               /**< originating separator of the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  feassetvars,        /**< variables in feasible set */
   int*                  nonfeassetvars,     /**< variables not in feasible set */
   int                   nfeassetvars,       /**< number of variables in feasible set */
   int                   nnonfeassetvars,    /**< number of variables not in feasible set */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   int* varsT1;
   int* varsT2;
   int* varsF;
   int* varsR;
   int* liftcoefs;
   SCIP_Real cutact;
   int nvarsT1;
   int nvarsT2;
   int nvarsF;
   int nvarsR;
   int liftrhs;
   int j;

   assert( cutoff != NULL );
   *cutoff = FALSE;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &varsT1, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsT2, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsF, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsR, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcoefs, nvars) );

   /* gets partition (T_1,T_2) of T, i.e. T_1 & T_2 = T and T_1 cap T_2 = emptyset, with T_1 not empty; chooses partition
    * as follows
    *   T_2 = { j in T : x*_j = 1 } and
    *   T_1 = T\T_2
    */
   getPartitionCovervars(scip, solvals, feassetvars, nfeassetvars, varsT1, varsT2, &nvarsT1, &nvarsT2);
   assert(nvarsT1 + nvarsT2 == nfeassetvars);

   /* changes partition (T_1,T_2) of feasible set T, if |T1| = 0, by moving one variable from T2 to T1 */
   if( nvarsT1 == 0 && nvarsT2 > 0)
   {
      SCIP_CALL( changePartitionFeasiblesetvars(scip, weights, varsT1, varsT2, &nvarsT1, &nvarsT2) );
      assert(nvarsT1 == 1);
   }
   assert(nvarsT2 == 0 || nvarsT1 > 0);

   /* gets partition (F,R) of N\T, i.e. F & R = N\T and F cap R = emptyset; chooses partition as follows
    *   R = { j in N\T : x*_j = 0 } and
    *   F = (N\T)\F
    */
   getPartitionNoncovervars(scip, solvals, nonfeassetvars, nnonfeassetvars, varsF, varsR, &nvarsF, &nvarsR);
   assert(nvarsF + nvarsR == nnonfeassetvars);
   assert(nvarsT1 + nvarsT2 + nvarsF + nvarsR == nvars - ntightened);

   /* sorts variables in F, T_2, and R according to the second level lifting sequence that will be used in the sequential
    * lifting procedure (the variable removed last from the initial cover does not have to be lifted first, therefore it
    * is included in the sorting routine)
    */
   SCIP_CALL( getLiftingSequence(scip, solvals, weights, varsF, varsT2, varsR, nvarsF, nvarsT2, nvarsR) );

   /* lifts extended weight inequality sum_{j in T_1} x_j <= |T_1| valid for
    *
    *    S^0 = { x in {0,1}^|T_1| : sum_{j in T_1} a_j x_j <= a_0 - sum_{j in T_2} a_j }
    *
    * to a valid inequality sum_{j in T_1} x_j + sum_{j in N\T_1} alpha_j x_j <= |T_1| + sum_{j in T_2} alpha_j for
    *
    *      S = { x in {0,1}^|N|   : sum_{j in N}   a_j x_j <= a_0 },
    *
    * uses sequential up-lifting for the variables in F, sequential down-lifting for the variable in T_2 and sequential
    * up-lifting for the variabels in R according to the second level lifting sequence
    */
   SCIP_CALL( sequentialUpAndDownLifting(scip, vars, nvars, ntightened, weights, capacity, solvals, varsT1, varsT2, varsF, varsR,
         nvarsT1, nvarsT2, nvarsF, nvarsR, nvarsT1, liftcoefs, &cutact, &liftrhs) );

   /* checks, if lifting yielded a violated cut */
   if( SCIPisEfficacious(scip, (cutact - liftrhs)/sqrt((SCIP_Real)MAX(liftrhs, 1))) )
   {
      SCIP_ROW* row;
      char name[SCIP_MAXSTRLEN];

      /* creates LP row */
      assert( cons == NULL || sepa == NULL );
      if( cons != NULL )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_ewseq%" SCIP_LONGINT_FORMAT "", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real)liftrhs,
               cons != NULL ? SCIPconsIsLocal(cons) : FALSE, FALSE,
               cons != NULL ? SCIPconsIsRemovable(cons) : TRUE) );
      }
      else if ( sepa != NULL )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_ewseq_%" SCIP_LONGINT_FORMAT "", SCIPsepaGetName(sepa), SCIPsepaGetNCutsFound(sepa));
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, FALSE, FALSE, TRUE) );
      }
      else
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nn_ewseq_%" SCIP_LONGINT_FORMAT "", *ncuts);
         SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, FALSE, FALSE, TRUE) );
      }

      /* adds all variables in the knapsack constraint with calculated lifting coefficient to the cut */
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      assert(nvarsT1 + nvarsT2 + nvarsF + nvarsR == nvars - ntightened);
      for( j = 0; j < nvarsT1; j++ )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsT1[j]], 1.0) );
      }
      for( j = 0; j < nvarsT2; j++ )
      {
         if( liftcoefs[varsT2[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsT2[j]], (SCIP_Real)liftcoefs[varsT2[j]]) );
         }
      }
      for( j = 0; j < nvarsF; j++ )
      {
         if( liftcoefs[varsF[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsF[j]], (SCIP_Real)liftcoefs[varsF[j]]) );
         }
      }
      for( j = 0; j < nvarsR; j++ )
      {
         if( liftcoefs[varsR[j]] > 0 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[varsR[j]], (SCIP_Real)liftcoefs[varsR[j]]) );
         }
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      /* checks, if cut is violated enough */
      if( SCIPisCutEfficacious(scip, sol, row) )
      {
         if( cons != NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
         (*ncuts)++;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &liftcoefs);
   SCIPfreeBufferArray(scip, &varsR);
   SCIPfreeBufferArray(scip, &varsF);
   SCIPfreeBufferArray(scip, &varsT2);
   SCIPfreeBufferArray(scip, &varsT1);

   return SCIP_OKAY;
}

/** separates lifted minimal cover inequalities using superadditive up-lifting for given knapsack problem */
static
SCIP_RETCODE separateSupLiftedMinimalCoverInequality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_SEPA*            sepa,               /**< originating separator of the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  mincovervars,       /**< mincover variables */
   int*                  nonmincovervars,    /**< nonmincover variables */
   int                   nmincovervars,      /**< number of mincover variables */
   int                   nnonmincovervars,   /**< number of nonmincover variables */
   SCIP_Longint          mincoverweight,     /**< weight of minimal cover */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_Real* realliftcoefs;
   SCIP_Real cutact;
   int liftrhs;

   assert( cutoff != NULL );
   *cutoff = FALSE;
   cutact = 0.0;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &realliftcoefs, nvars) );

   /* lifts minimal cover inequality sum_{j in C} x_j <= |C| - 1 valid for
    *
    *    S^0 = { x in {0,1}^|C| : sum_{j in C} a_j x_j <= a_0 }
    *
    * to a valid inequality sum_{j in C} x_j + sum_{j in N\C} alpha_j x_j <= |C| - 1 for
    *
    *      S = { x in {0,1}^|N|   : sum_{j in N}   a_j x_j <= a_0 },
    *
    * uses superadditive up-lifting for the variables in N\C.
    */
   SCIP_CALL( superadditiveUpLifting(scip, vars, nvars, ntightened, weights, capacity, solvals, mincovervars,
         nonmincovervars, nmincovervars, nnonmincovervars, mincoverweight, realliftcoefs, &cutact) );
   liftrhs = nmincovervars - 1;

   /* checks, if lifting yielded a violated cut */
   if( SCIPisEfficacious(scip, (cutact - liftrhs)/sqrt((SCIP_Real)MAX(liftrhs, 1))) )
   {
      SCIP_ROW* row;
      char name[SCIP_MAXSTRLEN];
      int j;

      /* creates LP row */
      assert( cons == NULL || sepa == NULL );
      if ( cons != NULL )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_mcsup%" SCIP_LONGINT_FORMAT "", SCIPconsGetName(cons), SCIPconshdlrGetNCutsFound(SCIPconsGetHdlr(cons)));
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real)liftrhs,
               cons != NULL ? SCIPconsIsLocal(cons) : FALSE, FALSE,
               cons != NULL ? SCIPconsIsRemovable(cons) : TRUE) );
      }
      else if ( sepa != NULL )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_mcsup%" SCIP_LONGINT_FORMAT "", SCIPsepaGetName(sepa), SCIPsepaGetNCutsFound(sepa));
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, FALSE, FALSE, TRUE) );
      }
      else
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nn_mcsup_%" SCIP_LONGINT_FORMAT "", *ncuts);
         SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)liftrhs, FALSE, FALSE, TRUE) );
      }

      /* adds all variables in the knapsack constraint with calculated lifting coefficient to the cut */
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      assert(nmincovervars + nnonmincovervars == nvars - ntightened);
      for( j = 0; j < nmincovervars; j++ )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[mincovervars[j]], 1.0) );
      }
      for( j = 0; j < nnonmincovervars; j++ )
      {
         assert(SCIPisFeasGE(scip, realliftcoefs[nonmincovervars[j]], 0.0));
         if( SCIPisFeasGT(scip, realliftcoefs[nonmincovervars[j]], 0.0) )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[nonmincovervars[j]], realliftcoefs[nonmincovervars[j]]) );
         }
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      /* checks, if cut is violated enough */
      if( SCIPisCutEfficacious(scip, sol, row) )
      {
         if( cons != NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
         (*ncuts)++;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &realliftcoefs);

   return SCIP_OKAY;
}

/** converts given cover C to a minimal cover by removing variables in the reverse order in which the variables were chosen
 *  to be in C, i.e. in the order of non-increasing (1 - x*_j)/a_j, if the transformed separation problem was used to find
 *  C and in the order of non-increasing (1 - x*_j), if the modified transformed separation problem was used to find C; 
 *  note that all variables with x*_j = 1 will be removed last
 */
static
SCIP_RETCODE makeCoverMinimal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< pointer to store cover variables */
   int*                  noncovervars,       /**< pointer to store noncover variables */
   int*                  ncovervars,         /**< pointer to store number of cover variables */
   int*                  nnoncovervars,      /**< pointer to store number of noncover variables */
   SCIP_Longint*         coverweight,        /**< pointer to store weight of cover */
   SCIP_Bool             modtransused        /**< TRUE if mod trans sepa prob was used to find cover */    
   )
{
   SORTKEYPAIR** sortkeypairs;
   SCIP_Longint minweight;
   int nsortkeypairs;
   int minweightidx;
   int j;
   int k;

   assert(scip != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars != NULL);
   assert(*ncovervars > 0);
   assert(nnoncovervars != NULL);
   assert(*nnoncovervars >= 0);
   assert(coverweight != NULL);
   assert(*coverweight > 0);
   assert(*coverweight > capacity);

   /* allocates temporary memory */
   nsortkeypairs = *ncovervars;
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeypairs, nsortkeypairs) );

   /* sorts C in the reverse order in which the variables were chosen to be in the cover, i.e. 
    *   such that (1 - x*_1)/a_1 >= ... >= (1 - x*_|C|)/a_|C|,  if          trans separation problem was used to find C 
    *   such that (1 - x*_1)     >= ... >= (1 - x*_|C|),        if modified trans separation problem was used to find C 
    * note that all variables with x*_j = 1 are in the end of the sorted C, so they will be removed last from C
    */
   assert(*ncovervars == nsortkeypairs);
   if( modtransused )
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         SCIP_CALL( SCIPallocBuffer(scip, &(sortkeypairs[j])) ); /*lint !e866 */

         sortkeypairs[j]->key1 = solvals[covervars[j]]; 
         sortkeypairs[j]->key2 = (SCIP_Real) weights[covervars[j]]; 
      }
   }
   else
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         SCIP_CALL( SCIPallocBuffer(scip, &(sortkeypairs[j])) ); /*lint !e866 */

         sortkeypairs[j]->key1 = (solvals[covervars[j]] - 1.0) / ((SCIP_Real) weights[covervars[j]]);
         sortkeypairs[j]->key2 = (SCIP_Real) (-weights[covervars[j]]);
      }
   }
   SCIPsortPtrInt((void**)sortkeypairs, covervars, compSortkeypairs, *ncovervars);

   /* gets j' with a_j' = min{ a_j : j in C } */
   minweightidx = 0;
   minweight = weights[covervars[minweightidx]];
   for( j = 1; j < *ncovervars; j++ )
   {
      if( weights[covervars[j]] <= minweight )
      {
         minweightidx = j;
         minweight = weights[covervars[minweightidx]];
      }
   }
   assert(minweightidx >= 0 && minweightidx < *ncovervars);
   assert(minweight > 0 && minweight <= *coverweight);

   j = 0;
   /* removes variables from C until the remaining variables form a minimal cover */
   while( j < *ncovervars && ((*coverweight) - minweight > capacity) )
   {
      assert(minweightidx >= j);
      assert(checkMinweightidx(weights, capacity, covervars, *ncovervars, *coverweight, minweightidx, j));

      /* if sum_{i in C} a_i - a_j <= a_0, j cannot be removed from C */
      if( (*coverweight) - weights[covervars[j]] <= capacity )
      {
	 ++j;
         continue;
      }

      /* adds j to N\C */
      noncovervars[*nnoncovervars] = covervars[j];
      (*nnoncovervars)++;

      /* removes j from C */
      (*coverweight) -= weights[covervars[j]];
      for( k = j; k < (*ncovervars) - 1; k++ )
         covervars[k] = covervars[k+1];
      (*ncovervars)--;

      /* updates j' with a_j' = min{ a_j : j in C } */
      if( j == minweightidx )
      {
         minweightidx = 0;
         minweight = weights[covervars[minweightidx]];
         for( k = 1; k < *ncovervars; k++ )
         {
            if( weights[covervars[k]] <= minweight )
            {
               minweightidx = k;
               minweight = weights[covervars[minweightidx]];
            }
         }
         assert(minweight > 0 && minweight <= *coverweight);
         assert(minweightidx >= 0 && minweightidx < *ncovervars);
      }
      else
      {
         assert(minweightidx > j);
         minweightidx--;
      }
      /* j needs to stay the same */
   }
   assert((*coverweight) > capacity);
   assert((*coverweight) - minweight <= capacity);

   /* frees temporary memory */
   for( j = nsortkeypairs-1; j >= 0; j-- )
      SCIPfreeBuffer(scip, &(sortkeypairs[j])); /*lint !e866 */
   SCIPfreeBufferArray(scip, &sortkeypairs);

   return SCIP_OKAY;
}

/** converts given initial cover C_init to a feasible set by removing variables in the reverse order in which
 *  they were chosen to be in C_init:
 *   non-increasing (1 - x*_j)/a_j,   if          transformed separation problem was used to find C_init
 *   non-increasing (1 - x*_j),       if modified transformed separation problem was used to find C_init.
 *  separates lifted extended weight inequalities using sequential up- and down-lifting for this feasible set
 *  and all subsequent feasible sets.
 */
static
SCIP_RETCODE getFeasibleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem */
   SCIP_SEPA*            sepa,               /**< originating separator of the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   int                   ntightened,         /**< number of variables with tightened upper bound */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< pointer to store cover variables */
   int*                  noncovervars,       /**< pointer to store noncover variables */
   int*                  ncovervars,         /**< pointer to store number of cover variables */
   int*                  nnoncovervars,      /**< pointer to store number of noncover variables */
   SCIP_Longint*         coverweight,        /**< pointer to store weight of cover */
   SCIP_Bool             modtransused,       /**< TRUE if mod trans sepa prob was used to find cover */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_Real* sortkeys;
   int j;
   int k;

   assert(scip != NULL);
   assert(covervars != NULL);
   assert(noncovervars != NULL);
   assert(ncovervars != NULL);
   assert(*ncovervars > 0);
   assert(nnoncovervars != NULL);
   assert(*nnoncovervars >= 0);
   assert(coverweight != NULL);
   assert(*coverweight > 0);
   assert(*coverweight > capacity);
   assert(*ncovervars + *nnoncovervars == nvars - ntightened);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortkeys, *ncovervars) );

   /* sorts C in the reverse order in which the variables were chosen to be in the cover, i.e.
    *   such that (1 - x*_1)/a_1 >= ... >= (1 - x*_|C|)/a_|C|,  if          trans separation problem was used to find C
    *   such that (1 - x*_1)     >= ... >= (1 - x*_|C|),        if modified trans separation problem was used to find C
    * note that all variables with x*_j = 1 are in the end of the sorted C, so they will be removed last from C
    */
   if( modtransused )
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         sortkeys[j] = solvals[covervars[j]];
         assert(SCIPisFeasGE(scip, sortkeys[j], 0.0));
      }
   }
   else
   {
      for( j = 0; j < *ncovervars; j++ )
      {
         sortkeys[j] = (solvals[covervars[j]] - 1.0) / ((SCIP_Real) weights[covervars[j]]);
         assert(SCIPisFeasLE(scip, sortkeys[j], 0.0));
      }
   }
   SCIPsortRealInt(sortkeys, covervars, *ncovervars);

   /* removes variables from C_init and separates lifted extended weight inequalities using sequential up- and down-lifting;
    * in addition to an extended weight inequality this gives cardinality inequalities */
   while( *ncovervars >= 2 )
   {
      /* adds first element of C_init to N\C_init */
      noncovervars[*nnoncovervars] = covervars[0];
      (*nnoncovervars)++;

      /* removes first element from C_init */
      (*coverweight) -= weights[covervars[0]];
      for( k = 0; k < (*ncovervars) - 1; k++ )
         covervars[k] = covervars[k+1];
      (*ncovervars)--;

      assert(*ncovervars + *nnoncovervars == nvars - ntightened);
      if( (*coverweight) <= capacity )
      {
         SCIP_CALL( separateSequLiftedExtendedWeightInequality(scip, cons, sepa, vars, nvars, ntightened, weights, capacity, solvals,
               covervars, noncovervars, *ncovervars, *nnoncovervars, sol, cutoff, ncuts) );
      }

      /* stop if cover is too large */
      if ( *ncovervars >= MAXCOVERSIZEITERLEWI )
         break;
   }

   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &sortkeys);

   return SCIP_OKAY;
}

/** separates different classes of valid inequalities for the 0-1 knapsack problem */
SCIP_RETCODE SCIPseparateKnapsackCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< originating constraint of the knapsack problem, or NULL */
   SCIP_SEPA*            sepa,               /**< originating separator of the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_SOL*             sol,                /**< primal SCIP solution to separate, NULL for current LP solution */
   SCIP_Bool             usegubs,            /**< should GUB information be used for separation? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_Real* solvals;
   int* covervars;
   int* noncovervars;
   SCIP_Bool coverfound;
   SCIP_Bool fractional;
   SCIP_Bool modtransused;
   SCIP_Longint coverweight;
   int ncovervars;
   int nnoncovervars;
   int ntightened;

   assert(scip != NULL);
   assert(capacity >= 0);
   assert(cutoff != NULL);
   assert(ncuts != NULL);

   *cutoff = FALSE;

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);
   assert(nvars > 0);
   assert(weights != NULL);

   /* increase age of constraint (age is reset to zero, if a cut was found) */
   if( cons != NULL )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* allocates temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &covervars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &noncovervars, nvars) );

   /* gets solution values of all problem variables */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, solvals) );

#ifdef SCIP_DEBUG
   {
      int i;

      SCIPdebugMsg(scip, "separate cuts for knapsack constraint originated by cons <%s>:\n",
         cons == NULL ? "-" : SCIPconsGetName(cons));
      for( i = 0; i < nvars; ++i )
      {
         SCIPdebugMsgPrint(scip, "%+" SCIP_LONGINT_FORMAT "<%s>(%g)", weights[i], SCIPvarGetName(vars[i]), solvals[i]);
      }
      SCIPdebugMsgPrint(scip, " <= %" SCIP_LONGINT_FORMAT "\n", capacity);
   }
#endif

   /* LMCI1 (lifted minimal cover inequalities using sequential up- and down-lifting) using GUB information
    */
   if( usegubs )
   {
      SCIP_GUBSET* gubset;

      SCIPdebugMsg(scip, "separate LMCI1-GUB cuts:\n");

      /* initializes partion of knapsack variables into nonoverlapping GUB constraints */
      SCIP_CALL( GUBsetCreate(scip, &gubset, nvars, weights, capacity) );

      /* constructs sophisticated partition of knapsack variables into nonoverlapping GUBs */
      SCIP_CALL( GUBsetGetCliquePartition(scip, gubset, vars, solvals) );
      assert(gubset->ngubconss <= nvars);

      /* gets a most violated initial cover C_init ( sum_{j in C_init} a_j > a_0 ) by using the
       * MODIFIED transformed separation problem and taking into account the following fixing:
       *   j in C_init,           if j in N_1 = {j in N : x*_j = 1} and
       *   j in N\C_init,         if j in N_0 = {j in N : x*_j = 0},
       * if one exists
       */
      modtransused = TRUE;
      SCIP_CALL( getCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars, &ncovervars,
            &nnoncovervars, &coverweight, &coverfound, modtransused, &ntightened, &fractional) );

      assert(!coverfound || !fractional || ncovervars + nnoncovervars == nvars - ntightened);

      /* if x* is not fractional we stop the separation routine */
      if( !fractional )
      {
         SCIPdebugMsg(scip, "   LMCI1-GUB terminated by no variable with fractional LP value.\n");

         /* frees memory for GUB set data structure */
         SCIP_CALL( GUBsetFree(scip, &gubset) );

         goto TERMINATE;
      }

      /* if no cover was found we stop the separation routine for lifted minimal cover inequality */
      if( coverfound )
      {
         /* converts initial cover C_init to a minimal cover C by removing variables in the reverse order in which the
          * variables were chosen to be in C_init; note that variables with x*_j = 1 will be removed last
          */
         SCIP_CALL( makeCoverMinimal(scip, weights, capacity, solvals, covervars, noncovervars, &ncovervars,
               &nnoncovervars, &coverweight, modtransused) );

         /* only separate with GUB information if we have at least one nontrivial GUB (with more than one variable) */
         if( gubset->ngubconss < nvars )
         {
            /* separates lifted minimal cover inequalities using sequential up- and down-lifting and GUB information */
            SCIP_CALL( separateSequLiftedMinimalCoverInequality(scip, cons, sepa, vars, nvars, ntightened, weights, capacity,
                  solvals, covervars, noncovervars, ncovervars, nnoncovervars, sol, gubset, cutoff, ncuts) );
         }
         else
         {
            /* separates lifted minimal cover inequalities using sequential up- and down-lifting, but do not use trivial
             * GUB information
             */
            SCIP_CALL( separateSequLiftedMinimalCoverInequality(scip, cons, sepa, vars, nvars, ntightened, weights, capacity,
                  solvals, covervars, noncovervars, ncovervars, nnoncovervars, sol, NULL, cutoff, ncuts) );
         }
      }

      /* frees memory for GUB set data structure */
      SCIP_CALL( GUBsetFree(scip, &gubset) );
   }
   else
   {
      /* LMCI1 (lifted minimal cover inequalities using sequential up- and down-lifting)
       * (and LMCI2 (lifted minimal cover inequalities using superadditive up-lifting))
       */

      /* gets a most violated initial cover C_init ( sum_{j in C_init} a_j > a_0 ) by using the
       * MODIFIED transformed separation problem and taking into account the following fixing:
       *   j in C_init,           if j in N_1 = {j in N : x*_j = 1} and
       *   j in N\C_init,         if j in N_0 = {j in N : x*_j = 0},
       * if one exists
       */
      SCIPdebugMsg(scip, "separate LMCI1 cuts:\n");
      modtransused = TRUE;
      SCIP_CALL( getCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars, &ncovervars,
            &nnoncovervars, &coverweight, &coverfound, modtransused, &ntightened, &fractional) );
      assert(!coverfound || !fractional || ncovervars + nnoncovervars == nvars - ntightened);

      /* if x* is not fractional we stop the separation routine */
      if( !fractional )
         goto TERMINATE;

      /* if no cover was found we stop the separation routine for lifted minimal cover inequality */
      if( coverfound )
      {
         /* converts initial cover C_init to a minimal cover C by removing variables in the reverse order in which the
          * variables were chosen to be in C_init; note that variables with x*_j = 1 will be removed last
          */
         SCIP_CALL( makeCoverMinimal(scip, weights, capacity, solvals, covervars, noncovervars, &ncovervars,
               &nnoncovervars, &coverweight, modtransused) );

         /* separates lifted minimal cover inequalities using sequential up- and down-lifting */
         SCIP_CALL( separateSequLiftedMinimalCoverInequality(scip, cons, sepa, vars, nvars, ntightened, weights, capacity,
               solvals, covervars, noncovervars, ncovervars, nnoncovervars, sol, NULL, cutoff, ncuts) );

         if( USESUPADDLIFT ) /*lint !e506 !e774*/
         {
            SCIPdebugMsg(scip, "separate LMCI2 cuts:\n");
            /* separates lifted minimal cover inequalities using superadditive up-lifting */
            SCIP_CALL( separateSupLiftedMinimalCoverInequality(scip, cons, sepa, vars, nvars, ntightened, weights, capacity,
                  solvals, covervars, noncovervars, ncovervars, nnoncovervars, coverweight, sol, cutoff, ncuts) );
         }
      }
   }

   /* LEWI (lifted extended weight inequalities using sequential up- and down-lifting) */
   if ( ! (*cutoff) )
   {
      /* gets a most violated initial cover C_init ( sum_{j in C_init} a_j > a_0 ) by using the
       * transformed separation problem and taking into account the following fixing:
       *   j in C_init,           if j in N_1 = {j in N : x*_j = 1} and
       *   j in N\C_init,         if j in N_0 = {j in N : x*_j = 0},
       * if one exists
       */
      SCIPdebugMsg(scip, "separate LEWI cuts:\n");
      modtransused = FALSE;
      SCIP_CALL( getCover(scip, vars, nvars, weights, capacity, solvals, covervars, noncovervars, &ncovervars,
            &nnoncovervars, &coverweight, &coverfound, modtransused, &ntightened, &fractional) );
      assert(fractional);
      assert(!coverfound || ncovervars + nnoncovervars == nvars - ntightened);

      /* if no cover was found we stop the separation routine */
      if( coverfound )
      {
         /* converts initial cover C_init to a feasible set by removing variables in the reverse order in which
          * they were chosen to be in C_init and separates lifted extended weight inequalities using sequential
          * up- and down-lifting for this feasible set and all subsequent feasible sets.
          */
         SCIP_CALL( getFeasibleSet(scip, cons, sepa, vars, nvars, ntightened, weights, capacity, solvals, covervars, noncovervars,
               &ncovervars, &nnoncovervars, &coverweight, modtransused, sol, cutoff, ncuts) );
      }
   }

 TERMINATE:
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &noncovervars);
   SCIPfreeBufferArray(scip, &covervars);
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}

/* relaxes given general linear constraint into a knapsack constraint and separates lifted knapsack cover inequalities */
SCIP_RETCODE SCIPseparateRelaxedKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< originating constraint of the knapsack problem, or NULL */
   SCIP_SEPA*            sepa,               /**< originating separator of the knapsack problem, or NULL */
   int                   nknapvars,          /**< number of variables in the continuous knapsack constraint */
   SCIP_VAR**            knapvars,           /**< variables in the continuous knapsack constraint */
   SCIP_Real*            knapvals,           /**< coefficients of the variables in the continuous knapsack constraint */
   SCIP_Real             valscale,           /**< -1.0 if lhs of row is used as rhs of c. k. constraint, +1.0 otherwise */
   SCIP_Real             rhs,                /**< right hand side of the continuous knapsack constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff was found */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR** consvars;
   SCIP_Real* binvals;
   SCIP_Longint* consvals;
   SCIP_Longint minact;
   SCIP_Longint maxact;
   SCIP_Real intscalar;
   SCIP_Bool success;
   int nbinvars;
   int nconsvars;
   int i;

   int* tmpindices;
   int tmp;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool noknapsackconshdlr;
   SCIP_Bool usegubs;

   assert(nknapvars > 0);
   assert(knapvars != NULL);
   assert(cutoff != NULL);

   tmpindices = NULL;

   SCIPdebugMsg(scip, "separate linear constraint <%s> relaxed to knapsack\n", cons != NULL ? SCIPconsGetName(cons) : "-");
   SCIPdebug( if( cons != NULL ) { SCIPdebugPrintCons(scip, cons, NULL); } );

   binvars = SCIPgetVars(scip);

   /* all variables which are of integral type can be potentially of binary type; this can be checked via the method SCIPvarIsBinary(var) */
   nbinvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   *cutoff = FALSE;

   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* set up data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );

   /* get conshdlrdata to use cleared memory */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      noknapsackconshdlr = TRUE;
      usegubs = DEFAULT_USEGUBS;

      SCIP_CALL( SCIPallocBufferArray(scip, &binvals, nbinvars) );
      BMSclearMemoryArray(binvals, nbinvars);
   }
   else
   {
      noknapsackconshdlr = FALSE;
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      usegubs = conshdlrdata->usegubs;

      SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices, nknapvars) );

      /* increase array size to avoid an endless loop in the next block; this might happen if continuous variables
       * change their types to SCIP_VARTYPE_BINARY during presolving
       */
      if( conshdlrdata->reals1size == 0 )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->reals1, conshdlrdata->reals1size, 1) );
         conshdlrdata->reals1size = 1;
         conshdlrdata->reals1[0] = 0.0;
      }

      assert(conshdlrdata->reals1size > 0);

      /* next if condition should normally not be true, because it means that presolving has created more binary
       * variables than binary + integer variables existed at the constraint initialization method, but for example if you would
       * transform all integers into their binary representation then it maybe happens
       */
      if( conshdlrdata->reals1size < nbinvars )
      {
         int oldsize = conshdlrdata->reals1size;

         conshdlrdata->reals1size = nbinvars;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->reals1, oldsize, conshdlrdata->reals1size) );
         BMSclearMemoryArray(&(conshdlrdata->reals1[oldsize]), conshdlrdata->reals1size - oldsize); /*lint !e866 */
      }
      binvals = conshdlrdata->reals1;

      /* check for cleared array, all entries have to be zero */
#ifndef NDEBUG
      for( tmp = nbinvars - 1; tmp >= 0; --tmp )
      {
         assert(binvals[tmp] == 0);
      }
#endif
   }

   tmp = 0;

   /* relax continuous knapsack constraint:
    * 1. make all variables binary:
    *    if x_j is continuous or integer variable substitute:
    *      - a_j < 0: x_j = lb  or  x_j = b*z + d with variable lower bound b*z + d with binary variable z
    *      - a_j > 0: x_j = ub  or  x_j = b*z + d with variable upper bound b*z + d with binary variable z
    * 2. convert coefficients of all variables to positive integers:
    *      - scale all coefficients a_j to a~_j integral
    *      - substitute  x~_j = 1 - x_j if a~_j < 0
    */

   /* replace integer and continuous variables with binary variables */
   for( i = 0; i < nknapvars; i++ )
   {
      SCIP_VAR* var;

      var = knapvars[i];

      if( SCIPvarIsBinary(var) && SCIPvarIsActive(var) )
      {
         SCIP_Real solval;
         assert(0 <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < nbinvars);

         solval = SCIPgetSolVal(scip, sol, var);

         /* knapsack relaxation assumes solution values between 0.0 and 1.0 for binary variables */
         if( SCIPisFeasLT(scip, solval, 0.0 )
               || SCIPisFeasGT(scip, solval, 1.0) )
         {
            SCIPdebugMsg(scip, "Solution value %.15g <%s> outside domain [0.0, 1.0]\n",
                  solval, SCIPvarGetName(var));
            goto TERMINATE;
         }

         binvals[SCIPvarGetProbindex(var)] += valscale * knapvals[i];
         if( !noknapsackconshdlr )
         {
            assert(tmpindices != NULL);

            tmpindices[tmp] = SCIPvarGetProbindex(var);
            ++tmp;
         }
         SCIPdebugMsg(scip, " -> binary variable %+.15g<%s>(%.15g)\n", valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var));
      }
      else if( valscale * knapvals[i] > 0.0 )
      {
         SCIP_VAR** zvlb;
         SCIP_Real* bvlb;
         SCIP_Real* dvlb;
         SCIP_Real bestlbsol;
         int bestlbtype;
         int nvlb;
         int j;

         /* a_j > 0: substitution with lb or vlb */
         nvlb = SCIPvarGetNVlbs(var);
         zvlb = SCIPvarGetVlbVars(var);
         bvlb = SCIPvarGetVlbCoefs(var);
         dvlb = SCIPvarGetVlbConstants(var);

         /* search for lb or vlb with maximal bound value */
         bestlbsol = SCIPvarGetLbGlobal(var);
         bestlbtype = -1;
         for( j = 0; j < nvlb; j++ )
         {
            /* use only numerical stable vlb with binary variable z */
            if( SCIPvarIsBinary(zvlb[j]) && SCIPvarIsActive(zvlb[j]) && REALABS(bvlb[j]) <= MAXABSVBCOEF )
            {
               SCIP_Real vlbsol;

               if( (bvlb[j] >= 0.0 && SCIPisGT(scip, bvlb[j] * SCIPvarGetLbLocal(zvlb[j]) + dvlb[j], SCIPvarGetUbLocal(var))) ||
                   (bvlb[j] <= 0.0 && SCIPisGT(scip, bvlb[j] * SCIPvarGetUbLocal(zvlb[j]) + dvlb[j], SCIPvarGetUbLocal(var))) )
               {
                  *cutoff = TRUE;
                  SCIPdebugMsg(scip, "variable bound <%s>[%g,%g] >= %g<%s>[%g,%g] + %g implies local cutoff\n",
                     SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                     bvlb[j], SCIPvarGetName(zvlb[j]), SCIPvarGetLbLocal(zvlb[j]), SCIPvarGetUbLocal(zvlb[j]), dvlb[j]);
                  goto TERMINATE;
               }

               assert(0 <= SCIPvarGetProbindex(zvlb[j]) && SCIPvarGetProbindex(zvlb[j]) < nbinvars);
               vlbsol = bvlb[j] * SCIPgetSolVal(scip, sol, zvlb[j]) + dvlb[j];
               if( SCIPisGE(scip, vlbsol, bestlbsol) )
               {
                  bestlbsol = vlbsol;
                  bestlbtype = j;
               }
            }
         }

         /* if no lb or vlb with binary variable was found, we have to abort */
         if( SCIPisInfinity(scip, -bestlbsol) )
            goto TERMINATE;

         if( bestlbtype == -1 )
         {
            rhs -= valscale * knapvals[i] * bestlbsol;
            SCIPdebugMsg(scip, " -> non-binary variable %+.15g<%s>(%.15g) replaced with lower bound %.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbGlobal(var), rhs);
         }
         else
         {
            assert(0 <= SCIPvarGetProbindex(zvlb[bestlbtype]) && SCIPvarGetProbindex(zvlb[bestlbtype]) < nbinvars);
            rhs -= valscale * knapvals[i] * dvlb[bestlbtype];
            binvals[SCIPvarGetProbindex(zvlb[bestlbtype])] += valscale * knapvals[i] * bvlb[bestlbtype];

            if( SCIPisInfinity(scip, REALABS(binvals[SCIPvarGetProbindex(zvlb[bestlbtype])])) )
               goto TERMINATE;

            if( !noknapsackconshdlr )
            {
               assert(tmpindices != NULL);

               tmpindices[tmp] = SCIPvarGetProbindex(zvlb[bestlbtype]);
               ++tmp;
            }
            SCIPdebugMsg(scip, " -> non-binary variable %+.15g<%s>(%.15g) replaced with variable lower bound %+.15g<%s>(%.15g) %+.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var),
               bvlb[bestlbtype], SCIPvarGetName(zvlb[bestlbtype]),
               SCIPgetSolVal(scip, sol, zvlb[bestlbtype]), dvlb[bestlbtype], rhs);
         }
      }
      else
      {
         SCIP_VAR** zvub;
         SCIP_Real* bvub;
         SCIP_Real* dvub;
         SCIP_Real bestubsol;
         int bestubtype;
         int nvub;
         int j;

         assert(valscale * knapvals[i] < 0.0);

         /* a_j < 0: substitution with ub or vub */
         nvub = SCIPvarGetNVubs(var);
         zvub = SCIPvarGetVubVars(var);
         bvub = SCIPvarGetVubCoefs(var);
         dvub = SCIPvarGetVubConstants(var);

         /* search for ub or vub with minimal bound value */
         bestubsol = SCIPvarGetUbGlobal(var);
         bestubtype = -1;
         for( j = 0; j < nvub; j++ )
         {
            /* use only numerical stable vub with active binary variable z */
            if( SCIPvarIsBinary(zvub[j]) && SCIPvarIsActive(zvub[j]) && REALABS(bvub[j]) <= MAXABSVBCOEF )
            {
               SCIP_Real vubsol;

               if( (bvub[j] >= 0.0 && SCIPisLT(scip, bvub[j] * SCIPvarGetUbLocal(zvub[j]) + dvub[j], SCIPvarGetLbLocal(var))) ||
                  (bvub[j] <= 0.0 && SCIPisLT(scip, bvub[j] * SCIPvarGetLbLocal(zvub[j]) + dvub[j], SCIPvarGetLbLocal(var))) )
               {
                  *cutoff = TRUE;
                  SCIPdebugMsg(scip, "variable bound <%s>[%g,%g] <= %g<%s>[%g,%g] + %g implies local cutoff\n",
                     SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                     bvub[j], SCIPvarGetName(zvub[j]), SCIPvarGetLbLocal(zvub[j]), SCIPvarGetUbLocal(zvub[j]), dvub[j]);
                  goto TERMINATE;
               }

               assert(0 <= SCIPvarGetProbindex(zvub[j]) && SCIPvarGetProbindex(zvub[j]) < nbinvars);
               vubsol = bvub[j] * SCIPgetSolVal(scip, sol, zvub[j]) + dvub[j];
               if( SCIPisLE(scip, vubsol, bestubsol) )
               {
                  bestubsol = vubsol;
                  bestubtype = j;
               }
            }
         }

         /* if no ub or vub with binary variable was found, we have to abort */
         if( SCIPisInfinity(scip, bestubsol) )
            goto TERMINATE;

         if( bestubtype == -1 )
         {
            rhs -= valscale * knapvals[i] * bestubsol;
            SCIPdebugMsg(scip, " -> non-binary variable %+.15g<%s>(%.15g) replaced with upper bound %.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetUbGlobal(var), rhs);
         }
         else
         {
            assert(0 <= SCIPvarGetProbindex(zvub[bestubtype]) && SCIPvarGetProbindex(zvub[bestubtype]) < nbinvars);
            rhs -= valscale * knapvals[i] * dvub[bestubtype];
            binvals[SCIPvarGetProbindex(zvub[bestubtype])] += valscale * knapvals[i] * bvub[bestubtype];

            if( SCIPisInfinity(scip, REALABS(binvals[SCIPvarGetProbindex(zvub[bestubtype])])) )
               goto TERMINATE;

            if( !noknapsackconshdlr )
            {
               assert(tmpindices != NULL);

               tmpindices[tmp] = SCIPvarGetProbindex(zvub[bestubtype]);
               ++tmp;
            }
            SCIPdebugMsg(scip, " -> non-binary variable %+.15g<%s>(%.15g) replaced with variable upper bound %+.15g<%s>(%.15g) %+.15g (rhs=%.15g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var),
               bvub[bestubtype], SCIPvarGetName(zvub[bestubtype]),
               SCIPgetSolVal(scip, sol, zvub[bestubtype]), dvub[bestubtype], rhs);
         }
      }
   }

   /* convert coefficients of all (now binary) variables to positive integers:
    *   - make all coefficients integral
    *   - make all coefficients positive (substitute negated variable)
    */
   nconsvars = 0;

   /* calculate scalar which makes all coefficients integral in relative allowed difference in between
    * -SCIPepsilon(scip) and KNAPSACKRELAX_MAXDELTA
    */
   SCIP_CALL( SCIPcalcIntegralScalar(binvals, nbinvars, -SCIPepsilon(scip), KNAPSACKRELAX_MAXDELTA,
         KNAPSACKRELAX_MAXDNOM, KNAPSACKRELAX_MAXSCALE, &intscalar, &success) );
   SCIPdebugMsg(scip, " -> intscalar = %.15g\n", intscalar);

   /* if coefficients cannot be made integral, we have to use a scalar of 1.0 and only round fractional coefficients down */
   if( !success )
      intscalar = 1.0;

   /* make all coefficients integral and positive:
    *  - scale a~_j = a_j * intscalar
    *  - substitute x~_j = 1 - x_j if a~_j < 0
    */
   rhs = rhs*intscalar;

   SCIPdebugMsg(scip, " -> rhs = %.15g\n", rhs);
   minact = 0;
   maxact = 0;
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_Longint val;

      val = (SCIP_Longint)SCIPfloor(scip, binvals[i]*intscalar);
      if( val == 0 )
         continue;

      if( val > 0 )
      {
         var = binvars[i];
         SCIPdebugMsg(scip, " -> positive scaled binary variable %+" SCIP_LONGINT_FORMAT "<%s> (unscaled %.15g): not changed (rhs=%.15g)\n",
            val, SCIPvarGetName(var), binvals[i], rhs);
      }
      else
      {
         assert(val < 0);

         SCIP_CALL( SCIPgetNegatedVar(scip, binvars[i], &var) );
         val = -val;
         rhs += val;
         SCIPdebugMsg(scip, " -> negative scaled binary variable %+" SCIP_LONGINT_FORMAT "<%s> (unscaled %.15g): substituted by (1 - <%s>) (rhs=%.15g)\n",
            -val, SCIPvarGetName(binvars[i]), binvals[i], SCIPvarGetName(var), rhs);
      }

      if( SCIPvarGetLbLocal(var) > 0.5 )
         minact += val;
      if( SCIPvarGetUbLocal(var) > 0.5 )
         maxact += val;
      consvals[nconsvars] = val;
      consvars[nconsvars] = var;
      nconsvars++;
   }

   if( nconsvars > 0 )
   {
      SCIP_Longint capacity;

      assert(consvars != NULL);
      assert(consvals != NULL);
      capacity = (SCIP_Longint)SCIPfeasFloor(scip, rhs);

#ifdef SCIP_DEBUG
      {
         SCIP_Real act;

         SCIPdebugMsg(scip, " -> linear constraint <%s> relaxed to knapsack:", cons != NULL ? SCIPconsGetName(cons) : "-");
         act = 0.0;
         for( i = 0; i < nconsvars; ++i )
         {
            SCIPdebugMsgPrint(scip, " %+" SCIP_LONGINT_FORMAT "<%s>(%.15g)", consvals[i], SCIPvarGetName(consvars[i]),
               SCIPgetSolVal(scip, sol, consvars[i]));
            act += consvals[i] * SCIPgetSolVal(scip, sol, consvars[i]);
         }
         SCIPdebugMsgPrint(scip, " <= %" SCIP_LONGINT_FORMAT " (%.15g) [act: %.15g, min: %" SCIP_LONGINT_FORMAT " max: %" SCIP_LONGINT_FORMAT "]\n",
            capacity, rhs, act, minact, maxact);
      }
#endif

      if( minact > capacity )
      {
         SCIPdebugMsg(scip, "minactivity of knapsack relaxation implies local cutoff\n");
         *cutoff = TRUE;
         goto TERMINATE;
      }

      if( maxact > capacity )
      {
         /* separate lifted cut from relaxed knapsack constraint */
         SCIP_CALL( SCIPseparateKnapsackCuts(scip, cons, sepa, consvars, nconsvars, consvals, capacity, sol, usegubs, cutoff, ncuts) );
      }
   }

 TERMINATE:
   /* free data structures */
   if( noknapsackconshdlr)
   {
      SCIPfreeBufferArray(scip, &binvals);
   }
   else
   {
      /* clear binvals */
      for( --tmp; tmp >= 0; --tmp)
      {
         assert(tmpindices != NULL);
         binvals[tmpindices[tmp]] = 0;
      }
      SCIPfreeBufferArray(scip, &tmpindices);
   }
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

/** separates given knapsack constraint */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_SOL*             sol,                /**< primal SCIP solution, NULL for current LP solution */
   SCIP_Bool             sepacuts,           /**< should knapsack cuts be separated? */
   SCIP_Bool             usegubs,            /**< should GUB information be used for separation? */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool violated;

   assert(ncuts != NULL);
   assert(cutoff != NULL);
   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "separating knapsack constraint <%s>\n", SCIPconsGetName(cons));

   /* check knapsack constraint itself for feasibility */
   SCIP_CALL( checkCons(scip, cons, sol, (sol != NULL), FALSE, &violated) );

   if( violated )
   {
      /* add knapsack constraint as LP row to the LP */
      SCIP_CALL( addRelaxation(scip, cons, cutoff) );
      (*ncuts)++;
   }
   else if( sepacuts )
   {
      SCIP_CALL( SCIPseparateKnapsackCuts(scip, cons, NULL, consdata->vars, consdata->nvars, consdata->weights,
            consdata->capacity, sol, usegubs, cutoff, ncuts) );
   }

   return SCIP_OKAY;
}

/** adds coefficient to constraint data */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_VAR*             var,                /**< variable to add to knapsack */
   SCIP_Longint          weight              /**< weight of variable in knapsack */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarIsBinary(var));
   assert(weight > 0);

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, (SCIP_Real)weight) );
   }

   /* check for fixed variable */
   if( SCIPvarGetLbGlobal(var) > 0.5 )
   {
      /* variable is fixed to one: reduce capacity */
      consdata->capacity -= weight;
   }
   else if( SCIPvarGetUbGlobal(var) > 0.5 )
   {
      SCIP_Bool negated;

      /* get binary representative of variable */
      SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &var, &negated) );

      /* insert coefficient */
      SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1, SCIPconsIsTransformed(cons)) );
      consdata->vars[consdata->nvars] = var;
      consdata->weights[consdata->nvars] = weight;
      consdata->nvars++;

      /* capture variable */
      SCIP_CALL( SCIPcaptureVar(scip, var) );

      /* install the rounding locks of variable */
      SCIP_CALL( lockRounding(scip, cons, var) );

      /* catch events */
      if( SCIPconsIsTransformed(cons) )
      {
         SCIP_CONSHDLRDATA* conshdlrdata;

         conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
         assert(conshdlrdata != NULL);
         SCIP_CALL( eventdataCreate(scip, &consdata->eventdata[consdata->nvars-1], cons, weight) );
         SCIP_CALL( SCIPcatchVarEvent(scip, var, EVENTTYPE_KNAPSACK,
               conshdlrdata->eventhdlr, consdata->eventdata[consdata->nvars-1],
               &consdata->eventdata[consdata->nvars-1]->filterpos) );

         if( !consdata->existmultaggr && SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_MULTAGGR )
            consdata->existmultaggr = TRUE;

         /* mark constraint to be propagated and presolved */
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
         consdata->presolvedtiming = 0;
         consdata->cliquesadded = FALSE; /* new coefficient might lead to larger cliques */
      }

      /* update weight sums */
      updateWeightSums(consdata, var, weight);

      consdata->sorted = FALSE;
      consdata->cliquepartitioned = FALSE;
      consdata->negcliquepartitioned = FALSE;
      consdata->merged = FALSE;
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* delete the coefficient from the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, -(SCIP_Real)consdata->weights[pos]) );
   }

   /* remove the rounding locks of variable */
   SCIP_CALL( unlockRounding(scip, cons, var) );

   /* drop events and mark constraint to be propagated and presolved */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
      assert(conshdlrdata != NULL);
      SCIP_CALL( SCIPdropVarEvent(scip, var, EVENTTYPE_KNAPSACK,
            conshdlrdata->eventhdlr, consdata->eventdata[pos], consdata->eventdata[pos]->filterpos) );
      SCIP_CALL( eventdataFree(scip, &consdata->eventdata[pos]) );

      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
      consdata->presolvedtiming = 0;
      consdata->sorted = (consdata->sorted && pos == consdata->nvars - 1);
   }

   /* decrease weight sums */
   updateWeightSums(consdata, var, -consdata->weights[pos]);

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->weights[pos] = consdata->weights[consdata->nvars-1];
   if( consdata->eventdata != NULL )
      consdata->eventdata[pos] = consdata->eventdata[consdata->nvars-1];

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   /* try to use old clique partitions */
   if( consdata->cliquepartitioned )
   {
      assert(consdata->cliquepartition != NULL);
      /* if the clique number is equal to the number of variables we have only cliques with one element, so we don't 
       * change the clique number */
      if( consdata->cliquepartition[consdata->nvars - 1] != consdata->nvars - 1 )
      {
         int oldcliqenum;

         oldcliqenum = consdata->cliquepartition[pos];
         consdata->cliquepartition[pos] = consdata->cliquepartition[consdata->nvars-1];

         /* the following if and else cases assure that we have increasing clique numbers */
         if( consdata->cliquepartition[pos] > pos )
            consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
         else
         {
            int i;
            int cliquenumbefore;

            /* if the old clique number was greater than the new one we have to check that before a bigger clique number 
             * occurs the same as the old one is still in the cliquepartition */
            if( oldcliqenum > consdata->cliquepartition[pos] )
            {
               for( i = 0; i < consdata->nvars; ++i )
                  if( oldcliqenum == consdata->cliquepartition[i] )
                     break;
                  else if( oldcliqenum < consdata->cliquepartition[i] )
                  {
                     consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
                     break;
                  }
               /* if we reached the end in the for loop, it means we have deleted the last element of the clique with
                * the biggest index, so decrease the number of cliques 
                */
               if( i == consdata->nvars )
                  --(consdata->ncliques);
            }
            /* if the old clique number was smaller than the new one we have to check the front for an element with 
             * clique number minus 1 */
            else if( oldcliqenum < consdata->cliquepartition[pos] )
            {
               cliquenumbefore = consdata->cliquepartition[pos] - 1;
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->cliquepartition[i] < cliquenumbefore; --i ); /*lint !e722*/

               if( i < cliquenumbefore )
                  consdata->cliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
            }
            /* if we deleted the last element of the clique with biggest index, we have to decrease the clique number */
            else if( pos == consdata->nvars - 1)
            {
               cliquenumbefore = consdata->cliquepartition[pos];
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->cliquepartition[i] < cliquenumbefore; --i ); /*lint !e722*/

               if( i < cliquenumbefore )
                  --(consdata->ncliques);
            }
            /* if the old clique number is equal to the new one the cliquepartition should be ok */
         }
      }
      else
         --(consdata->ncliques);
   }

   if( consdata->negcliquepartitioned )
   {
      assert(consdata->negcliquepartition != NULL);
      /* if the clique number is equal to the number of variables we have only cliques with one element, so we don't 
       * change the clique number */
      if( consdata->negcliquepartition[consdata->nvars-1] != consdata->nvars - 1 )
      {
         int oldcliqenum;

         oldcliqenum = consdata->negcliquepartition[pos];
         consdata->negcliquepartition[pos] = consdata->negcliquepartition[consdata->nvars-1];

         /* the following if and else cases assure that we have increasing clique numbers */
         if( consdata->negcliquepartition[pos] > pos )
            consdata->negcliquepartitioned = FALSE; /* recalculate the clique partition after a coefficient was removed */
         else
         {
            int i;
            int cliquenumbefore;

            /* if the old clique number was greater than the new one we have to check that, before a bigger clique number
             * occurs, the same as the old one occurs */
            if( oldcliqenum > consdata->negcliquepartition[pos] )
            {
               for( i = 0; i < consdata->nvars; ++i )
                  if( oldcliqenum == consdata->negcliquepartition[i] )
                     break;
                  else if( oldcliqenum < consdata->negcliquepartition[i] )
                  {
                     consdata->negcliquepartitioned = FALSE; /* recalculate the negated clique partition after a coefficient was removed */
                     break;
                  }
               /* if we reached the end in the for loop, it means we have deleted the last element of the clique with
                * the biggest index, so decrease the number of negated cliques 
                */
               if( i == consdata->nvars )
                  --(consdata->nnegcliques);
            }
            /* if the old clique number was smaller than the new one we have to check the front for an element with 
             * clique number minus 1 */
            else if( oldcliqenum < consdata->negcliquepartition[pos] )
            {
               cliquenumbefore = consdata->negcliquepartition[pos] - 1;
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->negcliquepartition[i] < cliquenumbefore; --i ); /*lint !e722*/

               if( i < cliquenumbefore )
                  consdata->negcliquepartitioned = FALSE; /* recalculate the negated clique partition after a coefficient was removed */
            }
            /* if we deleted the last element of the clique with biggest index, we have to decrease the clique number */
            else if( pos == consdata->nvars - 1)
            {
               cliquenumbefore = consdata->negcliquepartition[pos];
               for( i = pos - 1; i >= 0 && i >= cliquenumbefore && consdata->negcliquepartition[i] < cliquenumbefore; --i ); /*lint !e722*/ 

               if( i < cliquenumbefore )
                  --(consdata->nnegcliques);
            }
            /* otherwise if the old clique number is equal to the new one the cliquepartition should be ok */
         }
      }
      else
         --(consdata->nnegcliques);
   }

   --(consdata->nvars);

   return SCIP_OKAY;
}

/** removes all items with weight zero from knapsack constraint */
static
SCIP_RETCODE removeZeroWeights(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< knapsack constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( v = consdata->nvars-1; v >= 0; --v )
   {
      if( consdata->weights[v] == 0 )
      {
         SCIP_CALL( delCoefPos(scip, cons, v) );
      }
   }

   return SCIP_OKAY;
}

/* perform deletion of variables in all constraints of the constraint handler */
static
SCIP_RETCODE performVarDeletions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int v;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* iterate over all constraints */
   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      /* constraint is marked, that some of its variables were deleted */
      if( consdata->varsdeleted )
      {
         /* iterate over all variables of the constraint and delete them from the constraint */
         for( v = consdata->nvars - 1; v >= 0; --v )
         {
            if( SCIPvarIsDeleted(consdata->vars[v]) )
            {
               SCIP_CALL( delCoefPos(scip, conss[i], v) );
            }
         }
         consdata->varsdeleted = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** replaces multiple occurrences of a variable or its negation by a single coefficient */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   int v;
   int prev;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;

   if( consdata->merged ) 
      return SCIP_OKAY;

   if( consdata->nvars <= 1 )
   {
      consdata->merged = TRUE;
      return SCIP_OKAY;
   }

   assert(consdata->vars != NULL || consdata->nvars == 0);

   /* sorting array after indices of variables, that's only for faster merging */ 
   SCIPsortPtrPtrLongIntInt((void**)consdata->vars, (void**)consdata->eventdata, consdata->weights, 
      consdata->cliquepartition, consdata->negcliquepartition, SCIPvarCompActiveAndNegated, consdata->nvars);

   /* knapsack-sorting (decreasing weights) now lost */ 
   consdata->sorted = FALSE;

   v = consdata->nvars - 1;
   prev = v - 1;
   /* loop backwards through the items: deletion only affects rear items */
   while( prev >= 0 )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Bool negated1;
      SCIP_Bool negated2;

      negated1 = FALSE;
      negated2 = FALSE;

      var1 = consdata->vars[v];
      assert(SCIPvarIsBinary(var1));
      assert(SCIPvarIsActive(var1) || SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED);
      if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
      {
         var1 = SCIPvarGetNegatedVar(var1);
         negated1 = TRUE;
      }
      assert(var1 != NULL);

      var2 = consdata->vars[prev];
      assert(SCIPvarIsBinary(var2));
      assert(SCIPvarIsActive(var2) || SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED);
      if( SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED )
      {
         var2 = SCIPvarGetNegatedVar(var2);
         negated2 = TRUE;
      }
      assert(var2 != NULL);

      if( var1 == var2 )
      {
         /* both variables are either active or negated */
         if( negated1 == negated2 )
         {
            /* variables var1 and var2 are equal: add weight of var1 to var2, and delete var1 */
            consdataChgWeight(consdata, prev, consdata->weights[v] + consdata->weights[prev]);
            SCIP_CALL( delCoefPos(scip, cons, v) );
         }
         /* variables var1 and var2 are opposite: subtract smaller weight from larger weight, reduce capacity,
          * and delete item of smaller weight
          */
         else if( consdata->weights[v] == consdata->weights[prev] )
         {
            /* both variables eliminate themselves: w*x + w*(1-x) == w */
            consdata->capacity -= consdata->weights[v];
            SCIP_CALL( delCoefPos(scip, cons, v) ); /* this does not affect var2, because var2 stands before var1 */
            SCIP_CALL( delCoefPos(scip, cons, prev) );

            --prev;
         }
         else if( consdata->weights[v] < consdata->weights[prev] )
         {
            consdata->capacity -= consdata->weights[v];
            consdataChgWeight(consdata, prev, consdata->weights[prev] - consdata->weights[v]);
            assert(consdata->weights[prev] > 0);
            SCIP_CALL( delCoefPos(scip, cons, v) ); /* this does not affect var2, because var2 stands before var1 */
         }
         else
         {
            consdata->capacity -= consdata->weights[prev];
            consdataChgWeight(consdata, v, consdata->weights[v] - consdata->weights[prev]);
            assert(consdata->weights[v] > 0);
            SCIP_CALL( delCoefPos(scip, cons, prev) ); /* attention: normally we lose our order */
            /* restore order iff necessary */
            if( consdata->nvars != v ) /* otherwise the order still stands */
            {
               assert(prev == 0 || ((prev > 0) && (SCIPvarIsActive(consdata->vars[prev - 1]) || SCIPvarGetStatus(consdata->vars[prev - 1]) == SCIP_VARSTATUS_NEGATED)) );
               /* either that was the last pair or both, the negated and "normal" variable in front doesn't match var1, so the order is irrelevant */
               if( prev == 0 || (var1 != consdata->vars[prev - 1] && var1 != SCIPvarGetNegatedVar(consdata->vars[prev - 1])) )
                  --prev;
               else /* we need to let v at the same position*/
               {
                  consdata->cliquesadded = FALSE; /* reduced capacity might lead to larger cliques */
                  /* don't decrease v, the same variable may exist up front */
                  --prev;
                  continue;
               }
            }
         }
         consdata->cliquesadded = FALSE; /* reduced capacity might lead to larger cliques */
      }
      v = prev;
      --prev;
   }

   consdata->merged = TRUE;

   /* check infeasibility */
   if( consdata->onesweightsum > consdata->capacity )
   {
      SCIPdebugMsg(scip, "merge multiples detected cutoff.\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** in case the knapsack constraint is independent of every else, solve the knapsack problem (exactly) and apply the
 *  fixings (dual reductions)
 */
static
SCIP_RETCODE dualPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints  */
   SCIP_Bool*            deleted             /**< pointer to store if the constraint is deleted */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* profits;
   int* solitems;
   int* nonsolitems;
   int* items;
   SCIP_Real solval;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Bool applicable;
   int nsolitems;
   int nnonsolitems;
   int nvars;
   int v;

   assert(!SCIPconsIsModifiable(cons));

   /* constraints for which the check flag is set to FALSE, did not contribute to the lock numbers; therefore, we cannot
    * use the locks to decide for a dual reduction using this constraint; for example after a restart the cuts which are
    * added to the problems have the check flag set to FALSE
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   vars = consdata->vars;

   SCIP_CALL( SCIPallocBufferArray(scip, &profits, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &items, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, nvars) );

   applicable = TRUE;

   /* check if we can apply the dual reduction; this can be done if the knapsack has the only locks on this constraint;
    * collect object values which are the profits of the knapsack problem
    */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool negated;

      var = vars[v];
      assert(var != NULL);

      /* the variable should not be (globally) fixed */
      assert(SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5);

      if( SCIPvarGetNLocksDown(var) > 0 || SCIPvarGetNLocksUp(var) > 1 )
      {
         applicable = FALSE;
         break;
      }

      negated = FALSE;

      /* get the active variable */
      SCIP_CALL( SCIPvarGetProbvarBinary(&var, &negated) );
      assert(SCIPvarIsActive(var));

      if( negated )
         profits[v] = SCIPvarGetObj(var);
      else
         profits[v] = -SCIPvarGetObj(var);

      SCIPdebugMsg(scip, "variable <%s> -> item size %" SCIP_LONGINT_FORMAT ", profit <%g>\n",
         SCIPvarGetName(vars[v]), consdata->weights[v], profits[v]);
      items[v] = v;
   }

   if( applicable )
   {
      SCIP_Bool success;

      SCIPdebugMsg(scip, "the knapsack constraint <%s> is independent to rest of the problem\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);

      /* solve knapsack problem exactly */
      SCIP_CALL( SCIPsolveKnapsackExactly(scip, consdata->nvars, consdata->weights, profits, consdata->capacity,
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );

      if( success )
      {
         SCIP_VAR* var;

         /* apply solution of the knapsack as dual reductions */
         for( v = 0; v < nsolitems; ++v )
         {
            var = vars[solitems[v]];
            assert(var != NULL);

            SCIPdebugMsg(scip, "variable <%s> only locked up in knapsack constraints: dual presolve <%s>[%.15g,%.15g] >= 1.0\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
            SCIP_CALL( SCIPtightenVarLb(scip, var, 1.0, TRUE, &infeasible, &tightened) );
            assert(!infeasible);
            assert(tightened);
            (*nfixedvars)++;
         }

         for( v = 0; v < nnonsolitems; ++v )
         {
            var = vars[nonsolitems[v]];
            assert(var != NULL);

            SCIPdebugMsg(scip, "variable <%s> has no down locks: dual presolve <%s>[%.15g,%.15g] <= 0.0\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
            SCIP_CALL( SCIPtightenVarUb(scip, var, 0.0, TRUE, &infeasible, &tightened) );
            assert(!infeasible);
            assert(tightened);
            (*nfixedvars)++;
         }

         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         (*deleted) = TRUE;
      }
   }

   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &items);
   SCIPfreeBufferArray(scip, &profits);

   return SCIP_OKAY;
}

/** check if the knapsack constraint is parallel to objective function; if so update the cutoff bound and avoid that the
 *  constraint enters the LP by setting the initial and separated flag to FALSE
 */
static
SCIP_RETCODE checkParallelObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< knapsack constraint handler data */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real offset;
   SCIP_Real scale;
   SCIP_Real objval;
   SCIP_Bool applicable;
   SCIP_Bool negated;
   int nobjvars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   nobjvars = SCIPgetNObjVars(scip);

   /* check if the knapsack constraints has the same number of variables as the objective function and if the initial
    * and/or separated flag is set to FALSE
    */
   if( nvars != nobjvars || (!SCIPconsIsInitial(cons) && !SCIPconsIsSeparated(cons)) )
      return SCIP_OKAY;

   /* There are no variables in the ojective function and in the constraint. Thus, the constraint is redundant. Since we
    * have a pure feasibility problem, we do not want to set a cutoff or lower bound.
    */
   if( nobjvars == 0 )
      return SCIP_OKAY;

   vars = consdata->vars;
   assert(vars != NULL);

   applicable = TRUE;
   offset = 0.0;
   scale = 1.0;

   for( v = 0; v < nvars && applicable; ++v )
   {
      negated = FALSE;
      var = vars[v];
      assert(var != NULL);

      if( SCIPvarIsNegated(var) )
      {
         negated = TRUE;
         var = SCIPvarGetNegatedVar(var);
         assert(var != NULL);
      }

      objval = SCIPvarGetObj(var);

      /* if a variable has a zero objective coefficient the knapsack constraint is not parallel to objective function */
      if( SCIPisZero(scip, objval) )
         applicable = FALSE;
      else
      {
         SCIP_Real weight;

         weight = (SCIP_Real)consdata->weights[v];

         if( negated )
         {
            if( v == 0 )
            {
               /* the first variable defines the scale */
               scale = weight / -objval;

               offset += weight;
            }
            else if( SCIPisEQ(scip, -objval * scale, weight) )
               offset += weight;
            else
               applicable = FALSE;
         }
         else if( v == 0 )
         {
            /* the first variable define the scale */
            scale = weight / objval;
         }
         else if( !SCIPisEQ(scip, objval * scale, weight) )
            applicable = FALSE;
      }
   }

   if( applicable )
   {
      if( SCIPisPositive(scip, scale) && conshdlrdata->detectcutoffbound )
      {
         SCIP_Real cutoffbound;

         /* avoid that the knapsack constraint enters the LP since it is parallel to the objective function */
         SCIP_CALL( SCIPsetConsInitial(scip, cons, FALSE) );
         SCIP_CALL( SCIPsetConsSeparated(scip, cons, FALSE) );

         cutoffbound = (consdata->capacity - offset) / scale;

         /* increase the cutoff bound value by an epsilon to ensue that solution with the value of the cutoff bound are
          * still excepted
          */
         cutoffbound += SCIPcutoffbounddelta(scip);

         SCIPdebugMsg(scip, "constraint <%s> is parallel to objective function and provids a cutoff bound <%g>\n",
            SCIPconsGetName(cons), cutoffbound);

         if( cutoffbound < SCIPgetCutoffbound(scip) )
         {
            SCIPdebugMsg(scip, "update cutoff bound <%g>\n", cutoffbound);

            SCIP_CALL( SCIPupdateCutoffbound(scip, cutoffbound) );
         }
         else
         {
            /* in case the cutoff bound is worse then currently known one we avoid additionaly enforcement and
             * propagation
             */
            SCIP_CALL( SCIPsetConsEnforced(scip, cons, FALSE) );
            SCIP_CALL( SCIPsetConsPropagated(scip, cons, FALSE) );
         }
      }
      else if( SCIPisNegative(scip, scale) && conshdlrdata->detectlowerbound )
      {
         SCIP_Real lowerbound;

         /* avoid that the knapsack constraint enters the LP since it is parallel to the objective function */
         SCIP_CALL( SCIPsetConsInitial(scip, cons, FALSE) );
         SCIP_CALL( SCIPsetConsSeparated(scip, cons, FALSE) );

         lowerbound = (consdata->capacity - offset) / scale;

         SCIPdebugMsg(scip, "constraint <%s> is parallel to objective function and provids a lower bound <%g>\n",
            SCIPconsGetName(cons), lowerbound);

         SCIP_CALL( SCIPupdateLocalLowerbound(scip, lowerbound) );
      }
   }

   return SCIP_OKAY;
}

/** sort the variables and weights w.r.t. the clique partition; thereby ensure the current order of the variables when a
 *  weight of one variable is greater or equal another weight and both variables are in the same cliques */
static
SCIP_RETCODE stableSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< knapsack constraint data */
   SCIP_VAR**            vars,               /**< array for sorted variables */
   SCIP_Longint*         weights,            /**< array for sorted weights */
   int*                  cliquestartposs,    /**< starting position array for each clique */
   SCIP_Bool             usenegatedclique    /**< should negated or normal clique partition be used */
   )
{ 
   SCIP_VAR** origvars;
   int norigvars;
   SCIP_Longint* origweights;
   int* cliquepartition;
   int ncliques;

   SCIP_VAR*** varpointers;
   SCIP_Longint** weightpointers;
   int* cliquecount;

   int nextpos;
   int c;
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(vars != NULL);
   assert(weights != NULL);
   assert(cliquestartposs != NULL);

   origweights = consdata->weights;
   origvars = consdata->vars;
   norigvars = consdata->nvars;

   assert(origvars != NULL || norigvars == 0);
   assert(origweights != NULL || norigvars == 0);

   if( norigvars == 0 )
      return SCIP_OKAY;

   if( usenegatedclique )
   {
      assert(consdata->negcliquepartitioned);

      cliquepartition = consdata->negcliquepartition;
      ncliques = consdata->nnegcliques;
   }
   else
   {
      assert(consdata->cliquepartitioned);

      cliquepartition = consdata->cliquepartition;
      ncliques = consdata->ncliques;
   }

   assert(cliquepartition != NULL);
   assert(ncliques > 0);

   /* we first count all clique items and alloc temporary memory for a bucket sort */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquecount, ncliques) );
   BMSclearMemoryArray(cliquecount, ncliques);

   /* first we count for each clique the number of elements */
   for( v = norigvars - 1; v >= 0; --v )
   {
      assert(0 <= cliquepartition[v] && cliquepartition[v] < ncliques);
      ++(cliquecount[cliquepartition[v]]);
   }

   /*@todo: maybe it is better to put largest cliques up front */

#ifndef NDEBUG
   BMSclearMemoryArray(vars, norigvars);
   BMSclearMemoryArray(weights, norigvars);
#endif
   SCIP_CALL( SCIPallocBufferArray(scip, &varpointers, ncliques) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weightpointers, ncliques) );

   nextpos = 0;
   /* now we initialize all start pointers for each clique, so they will be ordered */
   for( c = 0; c < ncliques; ++c )
   {
      /* to reach the goal that all variables of each clique will be standing next to each other we will initialize the
       * starting pointers for each clique by adding the number of each clique to the last clique starting pointer
       * e.g. clique1 has 4 elements and clique2 has 3 elements the the starting pointer for clique1 will be the pointer
       *      to vars[0], the starting pointer to clique2 will be the pointer to vars[4] and to clique3 it will be 
       *      vars[7]
       *
       */
      varpointers[c] = (SCIP_VAR**) (vars + nextpos);
      cliquestartposs[c] = nextpos;
      weightpointers[c] = (SCIP_Longint*) (weights + nextpos);
      assert(cliquecount[c] > 0);
      nextpos += cliquecount[c];
      assert(nextpos > 0);
   }
   assert(nextpos == norigvars);
   cliquestartposs[c] = nextpos;

   /* now we copy all variable and weights to the right order */
   for( v = 0; v < norigvars; ++v )
   {
      *(varpointers[cliquepartition[v]]) = origvars[v];  /*lint !e613*/
      ++(varpointers[cliquepartition[v]]);
      *(weightpointers[cliquepartition[v]]) = origweights[v];  /*lint !e613*/
      ++(weightpointers[cliquepartition[v]]);
   }
#ifndef NDEBUG
   for( v = 0; v < norigvars; ++v )
   {
      assert(vars[v] != NULL);
      assert(weights[v] > 0);
   }
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &weightpointers);
   SCIPfreeBufferArray(scip, &varpointers);
   SCIPfreeBufferArray(scip, &cliquecount);

   return SCIP_OKAY;
}

/** deletes all fixed variables from knapsack constraint, and replaces variables with binary representatives */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off, or NULL if this
                                              *   information is not needed; in this case, we apply all fixings
                                              *   instead of stopping after the first infeasible one */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   if( cutoff != NULL )
      *cutoff = FALSE;

   SCIPdebugMsg(scip, "apply fixings:\n");
   SCIPdebugPrintCons(scip, cons, NULL);

   /* check infeasibility */
   if ( consdata->onesweightsum > consdata->capacity )
   {
      SCIPdebugMsg(scip, "apply fixings detected cutoff.\n");

      if( cutoff != NULL )
         *cutoff = TRUE;

      return SCIP_OKAY;
   }

   /* all multi-aggregations should be resolved */
   consdata->existmultaggr = FALSE;

   v = 0;
   while( v < consdata->nvars )
   {
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         consdata->capacity -= consdata->weights[v];
         SCIP_CALL( delCoefPos(scip, cons, v) );
         consdata->cliquesadded = FALSE; /* reduced capacity might lead to larger cliques */
      }
      else if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, v) );
      }
      else
      {
         SCIP_VAR* repvar;
         SCIP_VAR* negvar;
         SCIP_VAR* workvar;
         SCIP_Longint weight;
         SCIP_Bool negated;

         weight = consdata->weights[v];

         /* get binary representative of variable */
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );
         assert(repvar != NULL);

         /* check for multi-aggregation */
         if( SCIPvarIsNegated(repvar) )
         {
            workvar = SCIPvarGetNegatedVar(repvar);
            assert(workvar != NULL);
            negated = TRUE;
         }
         else
         {
            workvar = repvar;
            negated = FALSE;
         }

         /* @todo maybe resolve the problem that the eliminating of the multi-aggregation leads to a non-knapsack
          * constraint (converting into a linear constraint), for example the multi-aggregation consist of a non-binary
          * variable or due to resolving now their are non-integral coefficients or a non-integral capacity
          *
          * If repvar is not negated so workvar = repvar, otherwise workvar = 1 - repvar. This means,
          * weight * workvar = weight * (a_1*y_1 + ... + a_n*y_n + c)
          *
          * The explanation for  the following block:
          * 1a) If repvar is a multi-aggregated variable weight * repvar should be replaced by
          *     weight * (a_1*y_1 + ... + a_n*y_n + c).
          * 1b) If repvar is a negated variable of a multi-aggregated variable weight * repvar should be replaced by
          *     weight - weight * (a_1*y_1 + ... + a_n*y_n + c), for better further use here we switch the sign of weight
          *     so now we have the replacement -weight + weight * (a_1*y_1 + ... + a_n*y_n + c).
          * 2)  For all replacement variable we check:
          * 2a) weight * a_i < 0 than we add -weight * a_i * y_i_neg to the constraint and adjust the capacity through
          *     capacity -= weight * a_i caused by the negation of y_i.
          * 2b) weight * a_i >= 0 than we add weight * a_i * y_i to the constraint.
          * 3a) If repvar was not negated we need to subtract weight * c from capacity.
          * 3b) If repvar was negated we need to subtract weight * (c - 1) from capacity(note we switched the sign of
          *     weight in this case.
          */
         if( SCIPvarGetStatus(workvar) == SCIP_VARSTATUS_MULTAGGR )
         {
            SCIP_VAR** aggrvars;
            SCIP_Real* aggrscalars;
            SCIP_Real aggrconst;
            int naggrvars;
            int i;

            SCIP_CALL( SCIPflattenVarAggregationGraph(scip, workvar) );
            naggrvars = SCIPvarGetMultaggrNVars(workvar);
            aggrvars = SCIPvarGetMultaggrVars(workvar);
            aggrscalars = SCIPvarGetMultaggrScalars(workvar);
            aggrconst = SCIPvarGetMultaggrConstant(workvar);
            assert((aggrvars != NULL && aggrscalars != NULL) || naggrvars == 0);

            if( !SCIPisIntegral(scip, weight * aggrconst) )
            {
               SCIPerrorMessage("try to resolve a multi-aggregation with a non-integral value for weight*aggrconst = %g\n", weight*aggrconst);
               return SCIP_ERROR;
            }

            /* if workvar was negated, we have to flip the weight */
            if( negated )
               weight *= -1;

            for( i = naggrvars - 1; i >= 0; --i )
            {
               assert(aggrvars != NULL);
               assert(aggrscalars != NULL);

               if( !SCIPvarIsBinary(aggrvars[i]) )
               {
                  SCIPerrorMessage("try to resolve a multi-aggregation with a non-binary variable <%s>\n", aggrvars[i]);
                  return SCIP_ERROR;
               }
               if( !SCIPisIntegral(scip, weight * aggrscalars[i]) )
               {
                  SCIPerrorMessage("try to resolve a multi-aggregation with a non-integral value for weight*aggrscalars = %g\n", weight*aggrscalars[i]);
                  return SCIP_ERROR;
               }
               /* if the new coefficient is smaller than zero, we need to add the negated variable instead and adjust the capacity */
               if( SCIPisNegative(scip, weight * aggrscalars[i]) )
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, aggrvars[i], &negvar));
                  assert(negvar != NULL);
                  SCIP_CALL( addCoef(scip, cons, negvar, (SCIP_Longint)(SCIPfloor(scip, -weight * aggrscalars[i] + 0.5))) );
                  consdata->capacity -= (SCIP_Longint)(SCIPfloor(scip, weight * aggrscalars[i] + 0.5));
               }
               else
               {
                  SCIP_CALL( addCoef(scip, cons, aggrvars[i], (SCIP_Longint)(SCIPfloor(scip, weight * aggrscalars[i] + 0.5))) );
               }
            }
            /* delete old coefficient */
            SCIP_CALL( delCoefPos(scip, cons, v) );

            /* adjust the capacity with the aggregation constant and if necessary the extra weight through the negation */
            if( negated )
               consdata->capacity -= (SCIP_Longint)SCIPfloor(scip, weight * (aggrconst - 1) + 0.5);
            else
               consdata->capacity -= (SCIP_Longint)SCIPfloor(scip, weight * aggrconst + 0.5);

            if( consdata->capacity < 0 )
            {
               if( cutoff != NULL )
               {
                  *cutoff = TRUE;
                  break;
               }
            }
         }
         /* check, if the variable should be replaced with the representative */
         else if( repvar != var )
         {
            /* delete old (aggregated) variable */
            SCIP_CALL( delCoefPos(scip, cons, v) );

            /* add representative instead */
            SCIP_CALL( addCoef(scip, cons, repvar, weight) );
         }
         else
            ++v;
      }
   }
   assert(consdata->onesweightsum == 0);

   SCIPdebugMsg(scip, "after applyFixings, before merging:\n");
   SCIPdebugPrintCons(scip, cons, NULL);

   /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have to
    * clean up the constraint
    */
   if( cutoff != NULL && !(*cutoff) )
   {
      SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
      SCIPdebugMsg(scip, "after applyFixings and merging:\n");
      SCIPdebugPrintCons(scip, cons, NULL);
   }

   return SCIP_OKAY;
}


/** propagation method for knapsack constraints */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            redundant,          /**< pointer to store whether constraint is redundant */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   SCIP_Bool             usenegatedclique    /**< should negated clique information be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Longint* secondmaxweights;
   SCIP_Longint minweightsum;
   SCIP_Longint residualcapacity;

   int nvars;
   int i;
   int nnegcliques;

   SCIP_VAR** myvars;
   SCIP_Longint* myweights;
   int* cliquestartposs;
   int* cliqueendposs;
   SCIP_Longint localminweightsum;
   SCIP_Bool foundmax;
   int c;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(redundant != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;
   *redundant = FALSE;

   SCIPdebugMsg(scip, "propagating knapsack constraint <%s>\n", SCIPconsGetName(cons));

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

#ifndef NDEBUG
   /* assert that only active or negated variables are present */
   for( i = 0; i < consdata->nvars && consdata->merged; ++i )
   {
      assert(SCIPvarIsActive(consdata->vars[i]) || SCIPvarIsNegated(consdata->vars[i]) || SCIPvarGetStatus(consdata->vars[i]) == SCIP_VARSTATUS_FIXED);
   }
#endif

   usenegatedclique = usenegatedclique && consdata->merged;

   /* init for debugging */
   myvars = NULL;
   myweights = NULL;
   cliquestartposs = NULL;
   secondmaxweights = NULL;
   minweightsum = 0;
   nvars = consdata->nvars;
   /* make sure, the items are sorted by non-increasing weight */
   sortItems(consdata);

   do
   {
      localminweightsum = 0;

      /* (1) compute the minimum weight of the knapsack constraint using negated clique information;
       *     a negated clique means, that at most one of the clique variables can be zero
       *     - minweightsum = sum_{negated cliques C} ( sum(wi : i \in C) - W_max(C) ), where W_max(C) is the maximal weight of C
       *
       *     if for i \in C (a negated clique) oneweightsum + minweightsum - wi + W_max(C) > capacity => xi = 1
       *     since replacing i with the element of maximal weight leads to infeasibility
       */
      if( usenegatedclique && nvars > 0 )
      {
         SCIP_CONSHDLRDATA* conshdlrdata;
         conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
         assert(conshdlrdata != NULL);

         /* compute clique partitions */
         SCIP_CALL( calcCliquepartition(scip, conshdlrdata, consdata, FALSE, TRUE) );
         nnegcliques = consdata->nnegcliques;

         /* if we have no real negated cliques we can stop here */
         if( nnegcliques == nvars )
         {
            /* run the standard algorithm that does not involve cliques */
            usenegatedclique = FALSE;
            break;
         }

         /* allocate temporary memory and initialize it */
         SCIP_CALL( SCIPduplicateBufferArray(scip, &myvars, consdata->vars, nvars) );
         SCIP_CALL( SCIPduplicateBufferArray(scip, &myweights, consdata->weights, nvars) ) ;
         SCIP_CALL( SCIPallocBufferArray(scip, &cliquestartposs, nnegcliques + 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &cliqueendposs, nnegcliques) );
         SCIP_CALL( SCIPallocBufferArray(scip, &secondmaxweights, nnegcliques) );
         BMSclearMemoryArray(secondmaxweights, nnegcliques);

         /* resort variables to avoid quadratic algorithm later on */
         SCIP_CALL( stableSort(scip, consdata, myvars, myweights, cliquestartposs, TRUE) );

         /* save the end positions of the cliques because start positions are moved in the following loop */
         for( c = 0; c < nnegcliques; ++c )
         {
            cliqueendposs[c] = cliquestartposs[c+1] - 1;
            assert(cliqueendposs[c] - cliquestartposs[c] >= 0);
         }

         c = 0;
         foundmax = FALSE;
         i = 0;

         while( i < nvars )
         {
            /* ignore variables of the negated clique which are fixed to one since these are counted in
             * consdata->onesweightsum
             */

            /* if there are only one variable negated cliques left we can stop */
            if( nnegcliques - c == nvars - i )
            {
               minweightsum += localminweightsum;
               localminweightsum = 0;
               break;
            }

            /* for summing up the minimum active weights due to cliques we have to omit the biggest weights of each
             * clique, we can only skip this clique if this variables is not fixed to zero, otherwise we have to fix all
             * other clique variables to one
             */
            if( cliquestartposs[c] == i )
            {
               assert(myweights[i] > 0);
               ++c;
               minweightsum += localminweightsum;
               localminweightsum = 0;
               foundmax = TRUE;

               if( SCIPvarGetLbLocal(myvars[i]) > 0.5 )
                  foundmax = FALSE;

               if( SCIPvarGetUbLocal(myvars[i]) > 0.5 )
               {
                  ++i;
                  continue;
               }
            }

            if( SCIPvarGetLbLocal(myvars[i]) < 0.5 )
            {
               assert(myweights[i] > 0);

               if( SCIPvarGetUbLocal(myvars[i]) > 0.5 )
               {
                  assert(myweights[i] <= myweights[cliquestartposs[c - 1]]);

                  if( !foundmax )
                  {
                     foundmax = TRUE;

                     /* overwrite cliquestartpos to the position of the first unfixed variable in this clique */
                     cliquestartposs[c - 1] = i;
                     ++i;

                     continue;
                  }
                  /* memorize second max weight for each clique */
                  if( secondmaxweights[c - 1] == 0 )
                     secondmaxweights[c - 1] = myweights[i];

                  localminweightsum += myweights[i];
               }
               /* we found a fixed variable to zero so all other variables in this negated clique have to be fixed to one */
               else
               {
                  int v;
                  /* fix all other variables of the negated clique to 1 */
                  for( v = cliquestartposs[c - 1]; v < cliquestartposs[c]; ++v )
                  {
                     if( v != i && SCIPvarGetLbLocal(myvars[v]) < 0.5 )
                     {
                        SCIPdebugMsg(scip, " -> fixing variable <%s> to 1, due to negated clique information\n", SCIPvarGetName(myvars[v]));
                        SCIP_CALL( SCIPinferBinvarCons(scip, myvars[v], TRUE, cons, SCIPvarGetIndex(myvars[i]), &infeasible, &tightened) );

                        if( infeasible )
                        {
                           assert( SCIPvarGetUbLocal(myvars[v]) < 0.5 );

                           /* analyze the infeasibility if conflict analysis is applicable */
                           if( SCIPisConflictAnalysisApplicable(scip) )
                           {
                              /* conflict analysis can only be applied in solving stage */
                              assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPinProbing(scip));

                              /* initialize the conflict analysis */
                              SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

                              /* add the two variables which are fixed to zero within a negated clique */
                              SCIP_CALL( SCIPaddConflictBinvar(scip, myvars[i]) );
                              SCIP_CALL( SCIPaddConflictBinvar(scip, myvars[v]) );

                              /* start the conflict analysis */
                              SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
                           }
                           *cutoff = TRUE;
                           break;
                        }
                        assert(tightened);
                        ++(*nfixedvars);
                        SCIP_CALL( SCIPresetConsAge(scip, cons) );
                     }
                  }
                  if( *cutoff )
                     break;

                  /* reset local minweightsum for clique because all fixed to one variables are now counted in consdata->onesweightsum */
                  localminweightsum = 0;
                  /* we can jump to the end of this clique */
                  i = cliqueendposs[c - 1];
               }
            }
            ++i;
         }
         /* add last clique minweightsum */
         minweightsum += localminweightsum;

         SCIPdebugMsg(scip, "knapsack constraint <%s> has minimum weight sum of <%" SCIP_LONGINT_FORMAT ">\n",
            SCIPconsGetName(cons), minweightsum + consdata->onesweightsum );

         /* check, if weights of fixed variables don't exceeds knapsack capacity */
         if( !(*cutoff) && consdata->capacity >= minweightsum + consdata->onesweightsum )
         {
            SCIP_Longint maxcliqueweight = -1LL;

            /* loop over cliques */
            for( c = 0; c < nnegcliques; ++c )
            {
               SCIP_VAR* maxvar;
               SCIP_Bool maxvarfixed;
               int endvarposclique;
               int startvarposclique;

               assert(myvars != NULL);
               assert(nnegcliques == consdata->nnegcliques);
               assert(myweights != NULL);
               assert(secondmaxweights != NULL);
               assert(cliquestartposs != NULL);

               endvarposclique = cliqueendposs[c];
               startvarposclique = cliquestartposs[c];

               maxvar = myvars[startvarposclique];

               /* no need to process this negated clique because all variables are already fixed (which we detect from a fixed maxvar) */
               if( SCIPvarGetUbLocal(maxvar) - SCIPvarGetLbLocal(maxvar) < 0.5 )
                  continue;

               maxcliqueweight = myweights[startvarposclique];
               maxvarfixed = FALSE;
               /* if the sum of all weights of fixed variables to one plus the minimalweightsum (minimal weight which is already
                * used in this knapsack due to negated cliques) plus any weight minus the second largest weight in this clique
                * exceeds the capacity the maximum weight variable can be fixed to zero.
                */
               if( consdata->onesweightsum + minweightsum  + (maxcliqueweight - secondmaxweights[c]) > consdata->capacity )
               {
#ifndef NDEBUG
                  SCIP_Longint oldonesweightsum = consdata->onesweightsum;
#endif
                  assert(maxcliqueweight >= secondmaxweights[c]);
                  assert(SCIPvarGetLbLocal(maxvar) < 0.5 && SCIPvarGetUbLocal(maxvar) > 0.5);

                  SCIPdebugMsg(scip, " -> fixing variable <%s> to 0\n", SCIPvarGetName(maxvar));
                  SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  SCIP_CALL( SCIPinferBinvarCons(scip, maxvar, FALSE, cons, cliquestartposs[c], &infeasible, &tightened) );
                  assert(consdata->onesweightsum == oldonesweightsum);
                  assert(!infeasible);
                  assert(tightened);
                  (*nfixedvars)++;
                  maxvarfixed = TRUE;
               }
               /* the remaining cliques are singletons such that all subsequent variables have a weight that
                * fits into the knapsack
                */
               else if( nnegcliques - c == nvars - startvarposclique )
                  break;
               /* early termination of the remaining loop because no further variable fixings are possible:
                *
                * the gain in any of the following negated cliques (the additional term if the maximum weight variable was set to 1, and the second
                * largest was set to 0) does not suffice to infer additional variable fixings because
                *
                * - the cliques are sorted by decreasing maximum weight -> for all c' >= c: maxweights[c'] <= maxcliqueweight
                * - their second largest elements are at least as large as the smallest weight of the knapsack
                */
               else if( consdata->onesweightsum + minweightsum + (maxcliqueweight - consdata->weights[nvars - 1]) <= consdata->capacity )
                  break;

               /* loop over items with non-maximal weight (omitting the first position) */
               for( i = endvarposclique; i > startvarposclique; --i )
               {
                  /* there should be no variable fixed to 0 between startvarposclique + 1 and endvarposclique unless we
                   * messed up the clique preprocessing in the previous loop to filter those variables out */
                  assert(SCIPvarGetUbLocal(myvars[i]) > 0.5);

                  /* only check variables of negated cliques for which no variable is locally fixed */
                  if( SCIPvarGetLbLocal(myvars[i]) < 0.5 )
                  {
                     assert(maxcliqueweight >= myweights[i]);
                     assert(i == endvarposclique || myweights[i] >= myweights[i+1]);

                     /* we fix the members of this clique with non-maximal weight in two cases to 1:
                      *
                      * the maxvar was already fixed to 0 because it has a huge gain.
                      *
                      * if for i \in C (a negated clique) onesweightsum - wi + W_max(c) > capacity  => xi = 1
                      * since replacing i with the element of maximal weight leads to infeasibility */
                     if( maxvarfixed || consdata->onesweightsum + minweightsum - myweights[i] + maxcliqueweight > consdata->capacity  )
                     {
#ifndef NDEBUG
                        SCIP_Longint oldonesweightsum = consdata->onesweightsum;
#endif
                        SCIPdebugMsg(scip, " -> fixing variable <%s> to 1, due to negated clique information\n", SCIPvarGetName(myvars[i]));
                        SCIP_CALL( SCIPinferBinvarCons(scip, myvars[i], TRUE, cons, -i, &infeasible, &tightened) );
                        assert(consdata->onesweightsum == oldonesweightsum + myweights[i]);
                        assert(!infeasible);
                        assert(tightened);
                        ++(*nfixedvars);
                        SCIP_CALL( SCIPresetConsAge(scip, cons) );

                        /* update minweightsum because now the variable is fixed to one and its weight is counted by
                         * consdata->onesweightsum
                         */
                        minweightsum -= myweights[i];
                        assert(minweightsum >= 0);
                     }
                     else
                        break;
                  }
               }
#ifndef NDEBUG
               /* in debug mode, we assert that we did not miss possible fixings by the break above */
               for( ; i > startvarposclique; --i )
               {
                  SCIP_Bool varisfixed = SCIPvarGetUbLocal(myvars[i]) - SCIPvarGetLbLocal(myvars[i]) < 0.5;
                  SCIP_Bool exceedscapacity = consdata->onesweightsum + minweightsum - myweights[i] + maxcliqueweight > consdata->capacity;

                  assert(i == endvarposclique || myweights[i] >= myweights[i+1]);
                  assert(varisfixed || !exceedscapacity);
               }
#endif
            }
         }
         SCIPfreeBufferArray(scip, &secondmaxweights);
         SCIPfreeBufferArray(scip, &cliqueendposs);
         SCIPfreeBufferArray(scip, &cliquestartposs);
         SCIPfreeBufferArray(scip, &myweights);
         SCIPfreeBufferArray(scip, &myvars);
      }

      assert(consdata->negcliquepartitioned || minweightsum == 0);
   }
   while( FALSE );

   assert(usenegatedclique || minweightsum == 0);
   /* check, if weights of fixed variables already exceed knapsack capacity */
   if( consdata->capacity < minweightsum + consdata->onesweightsum )
   {
      SCIPdebugMsg(scip, " -> cutoff - fixed weight: %" SCIP_LONGINT_FORMAT ", capacity: %" SCIP_LONGINT_FORMAT " \n",
         consdata->onesweightsum, consdata->capacity);

      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;

      /* analyze the cutoff in SOLVING stage and if conflict analysis is turned on */
      if( (SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPinProbing(scip)) && SCIPisConflictAnalysisApplicable(scip) )
      {
         /* start conflict analysis with the fixed-to-one variables, add only as many as needed to exceed the capacity */
         SCIP_Longint weight;

         weight = 0;

         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

         for( i = 0; i < nvars && weight <= consdata->capacity; i++ )
         {
            if( SCIPvarGetLbLocal(consdata->vars[i]) > 0.5)
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
               weight += consdata->weights[i];
            }
         }

         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      return SCIP_OKAY;
   }

   /* the algorithm below is a special case of propagation involving negated cliques */
   if( !usenegatedclique )
   {
      assert(consdata->sorted);
      residualcapacity = consdata->capacity - consdata->onesweightsum;

      /* fix all variables to zero, that don't fit into the knapsack anymore */
      for( i = 0; i < nvars && consdata->weights[i] > residualcapacity; ++i )
      {
         /* if all weights of fixed variables to one plus any weight exceeds the capacity the variables have to be fixed
          * to zero
          */
         if( SCIPvarGetLbLocal(consdata->vars[i]) < 0.5 )
         {
            if( SCIPvarGetUbLocal(consdata->vars[i]) > 0.5 )
            {
               assert(consdata->onesweightsum + consdata->weights[i] > consdata->capacity);
               SCIPdebugMsg(scip, " -> fixing variable <%s> to 0\n", SCIPvarGetName(consdata->vars[i]));
               SCIP_CALL( SCIPresetConsAge(scip, cons) );
               SCIP_CALL( SCIPinferBinvarCons(scip, consdata->vars[i], FALSE, cons, i, &infeasible, &tightened) );
               assert(!infeasible);
               assert(tightened);
               (*nfixedvars)++;
            }
         }
      }
   }

   /* check if the knapsack is now redundant */
   if( !SCIPconsIsModifiable(cons) )
   {
      SCIP_Longint unfixedweightsum = consdata->onesweightsum;

      /* sum up the weights of all unfixed variables, plus the weight sum of all variables fixed to one already */
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarGetLbLocal(consdata->vars[i]) + 0.5 < SCIPvarGetUbLocal(consdata->vars[i]) )
         {
            unfixedweightsum += consdata->weights[i];

            /* the weight sum is larger than the capacity, so the constraint is not redundant */
            if( unfixedweightsum > consdata->capacity )
               return SCIP_OKAY;
         }
      }
      /* we summed up all (unfixed and fixed to one) weights and did not exceed the capacity, so the constraint is redundant */
      SCIPdebugMsg(scip, " -> knapsack constraint <%s> is redundant: weightsum=%" SCIP_LONGINT_FORMAT ", unfixedweightsum=%" SCIP_LONGINT_FORMAT ", capacity=%" SCIP_LONGINT_FORMAT "\n",
         SCIPconsGetName(cons), consdata->weightsum, unfixedweightsum, consdata->capacity);
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      *redundant = TRUE;
   }

   return SCIP_OKAY;
}

/** all but one variable fit into the knapsack constraint, so we can upgrade this constraint to an logicor constraint
 *  containing all negated variables of this knapsack constraint
 */
static
SCIP_RETCODE upgradeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  ndelconss,          /**< pointer to store the amount of deleted constraints */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONS* newcons;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 1);

   /* if the knapsack constraint consists only of two variables, we can upgrade it to a set-packing constraint */
   if( consdata->nvars == 2 )
   {
      SCIPdebugMsg(scip, "upgrading knapsack constraint <%s> to a set-packing constraint", SCIPconsGetName(cons));

      SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, SCIPconsGetName(cons), consdata->nvars, consdata->vars,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
   }
   /* if the knapsack constraint consists of at least three variables, we can upgrade it to a logicor constraint
    * containing all negated variables of the knapsack
    */
   else
   {
      SCIP_VAR** consvars;

      SCIPdebugMsg(scip, "upgrading knapsack constraint <%s> to a logicor constraint", SCIPconsGetName(cons));

      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, consdata->nvars) );
      SCIP_CALL( SCIPgetNegatedVars(scip, consdata->nvars, consdata->vars, consvars) );

      SCIP_CALL( SCIPcreateConsLogicor(scip, &newcons, SCIPconsGetName(cons), consdata->nvars, consvars,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );

      SCIPfreeBufferArray(scip, &consvars);
   }

   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   ++(*naddconss);

   SCIP_CALL( SCIPdelCons(scip, cons) );
   ++(*ndelconss);

   return SCIP_OKAY;
}

/** delete redundant variables
 *
 * i.e. 5x1 + 5x2 + 5x3 + 2x4 + 1x5 <= 13   =>   x4, x5 always fits into the knapsack, so we can delete them
 *
 * i.e. 5x1 + 5x2 + 5x3 + 2x4 + 1x5 <= 8 and we have the cliqueinformation (x1,x2,x3) is a clique
 *      =>   x4, x5 always fits into the knapsack, so we can delete them
 *
 * i.e. 5x1 + 5x2 + 5x3 + 1x4 + 1x5 <= 6 and we have the cliqueinformation (x1,x2,x3) is a clique and (x4,x5) too
 *      =>   we create the set partitioning constraint x4 + x5 <= 1 and delete them in this knapsack
 */
static
SCIP_RETCODE deleteRedundantVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_Longint          frontsum,           /**< sum of front items which fit if we try to take from the first till the last */
   int                   splitpos,           /**< split position till when all front items are fitting, splitpos is the
                                              *   first which did not fit */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the amount of changed sides */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Longint gcd;
   int nvars;
   int w;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(naddconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 < frontsum && frontsum < consdata->weightsum);
   assert(0 < splitpos && splitpos < consdata->nvars);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   vars = consdata->vars;
   weights = consdata->weights;
   nvars = consdata->nvars;
   capacity = consdata->capacity;

   /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
    * weight must not be sorted by their index
    */
#ifndef NDEBUG
   for( w = nvars - 1; w > 0; --w )
      assert(weights[w] <= weights[w-1]);
#endif

   /* if there are no variables rear to splitpos, the constraint has no redundant variables */
   if( consdata->nvars - 1 == splitpos )
      return SCIP_OKAY;

   assert(frontsum + weights[splitpos] > capacity);

   /* detect redundant variables */
   if( consdata->weightsum - weights[splitpos] <= capacity )
   {
      /* all rear items are redundant, because leaving one item in front and incl. of splitpos out the rear itmes always
       * fit
       */
      SCIPdebugMsg(scip, "Found redundant variables in constraint <%s>.\n", SCIPconsGetName(cons));

      /* delete items and update capacity */
      for( w = nvars - 1; w > splitpos; --w )
      {
         consdata->capacity -= weights[w];
         SCIP_CALL( delCoefPos(scip, cons, w) );
      }
      assert(w == splitpos);

      ++(*nchgsides);
      *nchgcoefs += (nvars - splitpos);

      /* division by greatest common divisor */
      gcd = weights[w];
      for( ; w >= 0 && gcd > 1; --w )
      {
         gcd = SCIPcalcGreComDiv(gcd, weights[w]);
      }

      /* normalize if possible */
      if( gcd > 1 )
      {
         for( w = splitpos; w >= 0; --w )
         {
            consdataChgWeight(consdata, w, weights[w]/gcd);
         }
         (*nchgcoefs) += nvars;

         consdata->capacity /= gcd;
         ++(*nchgsides);
      }

      /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
       * weight must not be sorted by their index
       */
#ifndef NDEBUG
      for( w = consdata->nvars - 1; w > 0; --w )
         assert(weights[w] <= weights[w - 1]);
#endif
   }
   /* rear items can only be redundant, when the sum is smaller to the weight at splitpos and all rear items would
    * always fit into the knapsack, therefor the item directly after splitpos needs to be smaller than the one at
    * splitpos and needs to fit into the knapsack
    */
   else if( conshdlrdata->disaggregation && frontsum + weights[splitpos + 1] <= capacity )
   {
      int* clqpart;
      int nclq;
      int len;

      len = nvars - (splitpos + 1);
      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &clqpart, len) );

      /* calculate clique partition */
      SCIP_CALL( SCIPcalcCliquePartition(scip, &(consdata->vars[splitpos+1]), len, clqpart, &nclq) );

      /* check if we found at least one clique */
      if( nclq < len )
      {
         SCIP_Longint maxactduetoclq;
         int cliquenum;

         maxactduetoclq = 0;
         cliquenum = 0;

         /* calculate maximum activity due to cliques */
         for( w = 0; w < len; ++w )
         {
            assert(clqpart[w] >= 0 && clqpart[w] <= w);
            if( clqpart[w] == cliquenum )
            {
               maxactduetoclq += weights[w + splitpos + 1];
               ++cliquenum;
            }
         }

         /* all rear items are redundant due to clique information, if maxactduetoclq is smaller than the weight before,
          * so delete them and create for all clique the corresponding clique constraints and update the capacity
          */
         if( frontsum + maxactduetoclq <= capacity )
         {
            SCIP_VAR** clqvars;
            int nclqvars;
            int c;

            assert(maxactduetoclq < weights[splitpos]);

            SCIPdebugMsg(scip, "Found redundant variables in constraint <%s> due to clique information.\n", SCIPconsGetName(cons));

            /* allocate temporary memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &clqvars, len - nclq + 1) );

            for( c = 0; c < nclq; ++c )
            {
               nclqvars = 0;
               for( w = 0; w < len; ++w )
               {
                  if( clqpart[w] == c )
                  {
                     clqvars[nclqvars] = vars[w + splitpos + 1];
                     ++nclqvars;
                  }
               }

               /* we found a real clique so extract this constraint, because we do not know who this information generated so */
               if( nclqvars > 1 )
               {
                  SCIP_CONS* cliquecons;
                  char name[SCIP_MAXSTRLEN];

                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clq_%" SCIP_LONGINT_FORMAT "_%d", SCIPconsGetName(cons), capacity, c);
                  SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, nclqvars, clqvars,
                        SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                        SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                        SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                        SCIPconsIsStickingAtNode(cons)) );
                  SCIPdebugMsg(scip, " -> adding clique constraint: ");
                  SCIPdebugPrintCons(scip, cliquecons, NULL);
                  SCIP_CALL( SCIPaddCons(scip, cliquecons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
                  ++(*naddconss);
               }
            }

            /* delete items and update capacity */
            for( w = nvars - 1; w > splitpos; --w )
            {
               SCIP_CALL( delCoefPos(scip, cons, w) );
               ++(*nchgcoefs);
            }
            consdata->capacity -= maxactduetoclq;
            assert(frontsum <= consdata->capacity);
            ++(*nchgsides);

            assert(w == splitpos);

            /* renew weights pointer */
            weights = consdata->weights;

            /* division by greatest common divisor */
            gcd = weights[w];
            for( ; w >= 0 && gcd > 1; --w )
            {
               gcd = SCIPcalcGreComDiv(gcd, weights[w]);
            }

            /* normalize if possible */
            if( gcd > 1 )
            {
               for( w = splitpos; w >= 0; --w )
               {
                  consdataChgWeight(consdata, w, weights[w]/gcd);
               }
               (*nchgcoefs) += nvars;

               consdata->capacity /= gcd;
               ++(*nchgsides);
            }

            /* free temporary memory */
            SCIPfreeBufferArray(scip, &clqvars);

            /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
             * weight must not be sorted by their index
             */
#ifndef NDEBUG
            for( w = consdata->nvars - 1; w > 0; --w )
               assert(weights[w] <= weights[w - 1]);
#endif
         }
      }

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &clqpart);
   }

   return SCIP_OKAY;
}

/* detect redundant variables which always fits into the knapsack
 *
 * i.e. 5x1 + 5x2 + 5x3 + 2x4 + 1x5 <= 13   =>   x4, x5 always fits into the knapsack, so we can delete them
 *
 * i.e. 5x1 + 5x2 + 5x3 + 2x4 + 1x5 <= 8 and we have the cliqueinformation (x1,x2,x3) is a clique
 *      =>   x4, x5 always fits into the knapsack, so we can delete them
 *
 * i.e. 5x1 + 5x2 + 5x3 + 1x4 + 1x5 <= 6 and we have the cliqueinformation (x1,x2,x3) is a clique and (x4,x5) too
 *      =>   we create the set partitioning constraint x4 + x5 <= 1 and delete them in this knapsack
 */
static
SCIP_RETCODE detectRedundantVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  ndelconss,          /**< pointer to store the amount of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the amount of changed sides */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Longint sum;
   int noldchgcoefs;
   int nvars;
   int v;
   int w;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(naddconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars >= 2);
   assert(consdata->weightsum > consdata->capacity);

   noldchgcoefs = *nchgcoefs;
   vars = consdata->vars;
   weights = consdata->weights;
   nvars = consdata->nvars;
   capacity = consdata->capacity;
   sum = 0;

   /* search for maximal fitting items */
   for( v = 0; v < nvars && sum + weights[v] <= capacity; ++v )
      sum += weights[v];

   assert(v < nvars);

   /* all but one variable fit into the knapsack, so we can upgrade this constraint to a logicor */
   if( v == nvars - 1 )
   {
      SCIP_CALL( upgradeCons(scip, cons, ndelconss, naddconss) );
      assert(SCIPconsIsDeleted(cons));

      return SCIP_OKAY;
   }

   if( v < nvars - 1 )
   {
      /* try to delete variables */
      SCIP_CALL( deleteRedundantVars(scip, cons, sum, v, nchgcoefs, nchgsides, naddconss) );
      assert(consdata->nvars > 1);

      /* all but one variable fit into the knapsack, so we can upgrade this constraint to a logicor */
      if( v == consdata->nvars - 1 )
      {
         SCIP_CALL( upgradeCons(scip, cons, ndelconss, naddconss) );
         assert(SCIPconsIsDeleted(cons));
      }

      return SCIP_OKAY;
   }

   /* if we already found some redundant variables, stop here */
   if( *nchgcoefs > noldchgcoefs )
      return SCIP_OKAY;

   assert(vars == consdata->vars);
   assert(weights == consdata->weights);
   assert(nvars == consdata->nvars);
   assert(capacity == consdata->capacity);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   /* calculate clique partition */
   SCIP_CALL( calcCliquepartition(scip, conshdlrdata, consdata, TRUE, FALSE) );

   /* check for real existing cliques */
   if( consdata->cliquepartition[v] < v )
   {
      SCIP_Longint sumfront;
      SCIP_Longint maxactduetoclqfront;
      int* clqpart;
      int cliquenum;


      sumfront = 0;
      maxactduetoclqfront = 0;

      clqpart = consdata->cliquepartition;
      cliquenum = 0;

      /* calculate maximal activity due to cliques */
      for( w = 0; w < nvars; ++w )
      {
         assert(clqpart[w] >= 0 && clqpart[w] <= w);
         if( clqpart[w] == cliquenum )
         {
            if( maxactduetoclqfront + weights[w] <= capacity )
            {
               maxactduetoclqfront += weights[w];
               ++cliquenum;
            }
            else
               break;
         }
         sumfront += weights[w];
      }
      assert(w >= v);

      /* if all items fit, then delete the whole constraint but create clique constraints which led to this
       * information
       */
      if( conshdlrdata->disaggregation && w == nvars )
      {
         SCIP_VAR** clqvars;
         int nclqvars;
         int c;
         int ncliques;

         assert(maxactduetoclqfront <= capacity);

         SCIPdebugMsg(scip, "Found redundant constraint <%s> due to clique information.\n", SCIPconsGetName(cons));

         ncliques = consdata->ncliques;

         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &clqvars, nvars - ncliques + 1) );

         for( c = 0; c < ncliques; ++c )
         {
            nclqvars = 0;
            for( w = 0; w < nvars; ++w )
            {
               if( clqpart[w] == c )
               {
                  clqvars[nclqvars] = vars[w];
                  ++nclqvars;
               }
            }

            /* we found a real clique so extract this constraint, because we do not know who this information generated so */
            if( nclqvars > 1 )
            {
               SCIP_CONS* cliquecons;
               char name[SCIP_MAXSTRLEN];

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clq_%" SCIP_LONGINT_FORMAT "_%d", SCIPconsGetName(cons), capacity, c);
               SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, nclqvars, clqvars,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                     SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                     SCIPconsIsStickingAtNode(cons)) );
               SCIPdebugMsg(scip, " -> adding clique constraint: ");
               SCIPdebugPrintCons(scip, cliquecons, NULL);
               SCIP_CALL( SCIPaddCons(scip, cliquecons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
               ++(*naddconss);
            }
         }

         /* delete old constraint */
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         ++(*ndelconss);

         SCIPfreeBufferArray(scip, &clqvars);

         return SCIP_OKAY;
      }

      if( w > v && w < nvars - 1 )
      {
         /* try to delete variables */
         SCIP_CALL( deleteRedundantVars(scip, cons, sumfront, w, nchgcoefs, nchgsides, naddconss) );
      }
   }

   return SCIP_OKAY;
}

/** divides weights by their greatest common divisor and divides capacity by the same value, rounding down the result */
static
void normalizeWeights(
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Longint gcd;
   int i;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL); /* we are in presolve, so no LP row exists */
   assert(consdata->onesweightsum == 0); /* all fixed variables should have been removed */
   assert(consdata->weightsum > consdata->capacity); /* otherwise, the constraint is redundant */
   assert(consdata->nvars >= 1);

   /* sort items, because we can stop earlier if the smaller weights are evaluated first */
   sortItems(consdata);

   gcd = consdata->weights[consdata->nvars-1];
   for( i = consdata->nvars-2; i >= 0 && gcd >= 2; --i )
   {
      assert(SCIPvarGetLbLocal(consdata->vars[i]) < 0.5);
      assert(SCIPvarGetUbLocal(consdata->vars[i]) > 0.5); /* all fixed variables should have been removed */

      gcd = SCIPcalcGreComDiv(gcd, consdata->weights[i]);
   }

   if( gcd >= 2 )
   {
      SCIPdebugMessage("knapsack constraint <%s>: dividing weights by %" SCIP_LONGINT_FORMAT "\n", SCIPconsGetName(cons), gcd);

      for( i = 0; i < consdata->nvars; ++i )
      {
         consdataChgWeight(consdata, i, consdata->weights[i]/gcd);
      }
      consdata->capacity /= gcd;
      (*nchgcoefs) += consdata->nvars;
      (*nchgsides)++;

      /* weight should still be sorted, because the reduction preserves this */
#ifndef NDEBUG
      for( i = consdata->nvars - 1; i > 0; --i )
         assert(consdata->weights[i] <= consdata->weights[i - 1]);
#endif
      consdata->sorted = TRUE;
   }
}

/** dual weights tightening for knapsack constraints
 *
 *  1. a) check if all two pairs exceed the capacity, then we can upgrade this constraint to a set-packing constraint
 *     b) check if all but the smallest weight fit into the knapsack,  then we can upgrade this constraint to a logicor
 *        constraint
 *
 *  2. check if besides big coefficients, that fit only by itself, for a certain amount of variables all combination of
 *     these are a minimal cover, then might reduce the weights and the capacity, e.g.
 *
 *     +219y1 + 180y2 + 74x1 + 70x2 + 63x3 + 62x4 + 53x5 <= 219  <=>  3y1 + 3y2 + x1 + x2 + x3 + x4 + x5 <= 3
 *
 *  3. use the duality between a^Tx <= capacity   <=>   a^T~x >= weightsum - capacity to tighten weights, e.g.
 *
 *     11x1 + 10x2 + 7x3 + 7x4 + 5x5 <= 27    <=>   11~x1 + 10~x2 + 7~x3 + 7~x4 + 5~x5 >= 13
 *
 *     the above constraint can be changed to       8~x1 + 8~x2 + 6.5~x3 + 6.5~x4 + 5~x5 >= 13
 *
 *     16~x1 + 16~x2 + 13~x3 + 13~x4 + 10~x5 >= 26   <=>   16x1 + 16x2 + 13x3 + 13x4 + 10x5 <= 42
 */
static
SCIP_RETCODE dualWeightsTightening(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  ndelconss,          /**< pointer to store the amount of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the amount of changed sides */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Longint* weights;
   SCIP_Longint dualcapacity;
   SCIP_Longint reductionsum;
   SCIP_Longint capacity;
   SCIP_Longint exceedsum;
   int oldnchgcoefs;
   int nvars;
   int vbig;
   int v;
   int w;
#ifndef NDEBUG
   int oldnchgsides;
#endif

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(naddconss != NULL);

#ifndef NDEBUG
   oldnchgsides = *nchgsides;
#endif

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->weightsum > consdata->capacity);
   assert(consdata->nvars >= 2);
   assert(consdata->sorted);

   /* constraint should be merged */
   assert(consdata->merged);

   nvars = consdata->nvars;
   weights = consdata->weights;
   capacity = consdata->capacity;

   oldnchgcoefs = *nchgcoefs;

   /* case 1. */
   if( weights[nvars - 1] + weights[nvars - 2] > capacity )
   {
      SCIP_CONS* newcons;

      /* two variable are enough to exceed the constraint, so we can update it to a set-packing
       *
       * e.g. 5x1 + 4x2 + 3x3 <= 5   <=>    x1 + x2 + x3 <= 1
       */
      SCIPdebugMsg(scip, "upgrading knapsack constraint <%s> to a set-packing constraint", SCIPconsGetName(cons));

      SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, SCIPconsGetName(cons), consdata->nvars, consdata->vars,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );

      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      ++(*naddconss);

      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);

      return SCIP_OKAY;
   }

   /* all but one variable fit into the knapsack, so we can upgrade this constraint to a logicor */
   if( consdata->weightsum - weights[nvars - 1] <= consdata->capacity )
   {
      SCIP_CALL( upgradeCons(scip, cons, ndelconss, naddconss) );
      assert(SCIPconsIsDeleted(cons));

      return SCIP_OKAY;
   }

   /* early termination, if the pair with biggest coeffcients together does not exceed the dualcapacity */
   /* @todo might be changed/removed when improving the coeffcients tightening */
   if( consdata->weightsum - capacity > weights[0] + weights[1] )
      return SCIP_OKAY;

   /* case 2. */

   v = 0;

   /* @todo generalize the following algorithm for several parts of the knapsack
    *
    * the following is done without looking at the dualcapacity; it is enough to check whether for a certain amount of
    * variables each combination is a minimal cover, some examples
    *
    * +74x1 + 70x2 + 63x3 + 62x4 + 53x5 <= 219     <=>    74~x1 + 70~x2 + 63~x3 + 62~x4 + 53~x5 >= 103
    *                                              <=>      ~x1 +   ~x2 +   ~x3 +   ~x4 +   ~x5 >= 2
    *                                              <=>       x1 +    x2 +    x3 +    x4 +    x5 <= 3
    *
    * +219y1 + 180y_2 +74x1 + 70x2 + 63x3 + 62x4 + 53x5 <= 219  <=>  3y1 + 3y2 + x1 + x2 + x3 + x4 + x5 <= 3
    *
    */

   /* determine big weights that fit only by itself */
   while( v < nvars && weights[v] + weights[nvars - 1] > capacity )
      ++v;

   vbig = v;
   assert(vbig < nvars - 1);
   exceedsum = 0;

   /* determine the amount needed to exceed the capacity */
   while( v < nvars && exceedsum <= capacity )
   {
      exceedsum += weights[v];
      ++v;
   }

   /* if we exceeded the capacity we might reduce the weights */
   if( exceedsum > capacity )
   {
      assert(vbig > 0 || v < nvars);

      /* all small weights were needed to exceed the capacity */
      if( v == nvars )
      {
         SCIP_Longint newweight = (SCIP_Longint)nvars - vbig - 1;
         assert(newweight > 0);

         /* reduce big weights */
         for( v = 0; v < vbig; ++v )
         {
            if( weights[v] > newweight )
            {
               consdataChgWeight(consdata, v, newweight);
               ++(*nchgcoefs);
            }
         }

         /* reduce small weights */
         for( ; v < nvars; ++v )
         {
            if( weights[v] > 1 )
            {
               consdataChgWeight(consdata, v, 1LL);
               ++(*nchgcoefs);
            }
         }

         consdata->capacity = newweight;

         /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
          * weight must not be sorted by their index
          */
#ifndef NDEBUG
         for( v = nvars - 1; v > 0; --v )
            assert(weights[v] <= weights[v-1]);
#endif

         return SCIP_OKAY;
      }
      /* a certain amount of small variables exceed the capacity, so check if this holds for all combinations of the
       * small weights
       */
      else
      {
         SCIP_Longint exceedsumback = 0;
         int nexceed = v - vbig;

         assert(nexceed > 1);

         /* determine weightsum of the same amount as before but of the smallest weight */
         for( w = nvars - 1; w >= nvars - nexceed; --w )
            exceedsumback += weights[w];

         assert(w >= 0);

         /* if the same amount but with the smallest possible weights also exceed the capacity, it holds for all
          * combinations of all small weights
          */
         if( exceedsumback > capacity )
         {
            SCIP_Longint newweight = nexceed - 1;

            /* taking out the smallest element needs to fit */
            assert(exceedsumback - weights[nvars - 1] <= capacity);

            /* reduce big weights */
            for( v = 0; v < vbig; ++v )
            {
               if( weights[v] > newweight )
               {
                  consdataChgWeight(consdata, v, newweight);
                  ++(*nchgcoefs);
               }
            }

            /* reduce small weights */
            for( ; v < nvars; ++v )
            {
               if( weights[v] > 1 )
               {
                  consdataChgWeight(consdata, v, 1LL);
                  ++(*nchgcoefs);
               }
            }

            consdata->capacity = newweight;

            /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
             * weight must not be sorted by their index
             */
#ifndef NDEBUG
            for( v = nvars - 1; v > 0; --v )
               assert(weights[v] <= weights[v-1]);
#endif
            return SCIP_OKAY;
         }
      }
   }
   else
   {
      /* if the following assert fails we have either a redundant constraint or a set-packing constraint, this should
       * not happen here
       */
      assert(vbig > 0 && vbig < nvars);

      /* either choose a big coefficients or all other variables
       *
       * 973x1 + 189x2 + 189x3 + 145x4 + 110x5 + 104x6 + 93x7 + 71x8 + 68x9 + 10x10 <= 979
       *
       * either choose x1, or all other variables (weightsum of x2 to x10 is 979 above), so we can tighten this
       * constraint to
       *
       * 9x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 <= 9
       */

      if( weights[vbig - 1] > (SCIP_Longint)nvars - vbig || weights[vbig] > 1 )
      {
         SCIP_Longint newweight = (SCIP_Longint)nvars - vbig;
#ifndef NDEBUG
         SCIP_Longint resweightsum = consdata->weightsum;

         for( v = 0; v < vbig; ++v )
            resweightsum -= weights[v];

         assert(exceedsum == resweightsum);
#endif
         assert(newweight > 0);

         /* reduce big weights */
         for( v = 0; v < vbig; ++v )
         {
            if( weights[v] > newweight )
            {
               consdataChgWeight(consdata, v, newweight);
               ++(*nchgcoefs);
            }
         }

         /* reduce small weights */
         for( ; v < nvars; ++v )
         {
            if( weights[v] > 1 )
            {
               consdataChgWeight(consdata, v, 1LL);
               ++(*nchgcoefs);
            }
         }

         consdata->capacity = newweight;

         /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
          * weight must not be sorted by their index
          */
#ifndef NDEBUG
         for( v = nvars - 1; v > 0; --v )
            assert(weights[v] <= weights[v-1]);
#endif
         return SCIP_OKAY;
      }
   }

   /* case 3. */

   dualcapacity = consdata->weightsum - capacity;
   reductionsum = 0;
   v = 0;

   /* reduce big weights
    *
    * e.g. 11x0 + 11x1 + 10x2 + 10x3 <= 32   <=>    11~x0 + 11~x1 + 10~x2 + 10~x3 >= 10
    *                                        <=>    10~x0 + 10~x1 + 10~x2 + 10~x3 >= 10
    *                                        <=>       x0 +    x1 +    x2 +    x3 <= 3
    */
   while( weights[v] > dualcapacity )
   {
      reductionsum += (weights[v] - dualcapacity);
      consdataChgWeight(consdata, v, dualcapacity);
      ++v;
      assert(v < nvars);
   }
   (*nchgcoefs) += v;

   /* skip weights equal to the dualcapacity, because we cannot change them  */
   while( v < nvars && weights[v] == dualcapacity )
      ++v;

   /* any negated variable out of the first n - 1 items is enough to fulfill the constraint, so we can update it to a logicor
    * after a possible removal of the last, redundant item
    *
    * e.g. 10x1 + 10x2 + 10x3 <= 20   <=>    10~x1 + 10~x2 + 10~x3 >= 10  <=>   ~x1 + ~x2 + ~x3 >= 1
    */
   if( v >= nvars - 1 )
   {
      /* the last weight is not enough to satisfy the dual capacity -> remove this redundant item */
      if( v == nvars - 1 )
      {
         SCIP_CALL( delCoefPos(scip, cons, nvars - 1) );
      }
      SCIP_CALL( upgradeCons(scip, cons, ndelconss, naddconss) );
      assert(SCIPconsIsDeleted(cons));

      return SCIP_OKAY;
   }
   else /* v < nvars - 1 <=> at least two items with weight smaller than the dual capacity */
   {
      /* @todo generalize the following algorithm for more than two variables */

      if( weights[nvars - 1] + weights[nvars - 2] >= dualcapacity )
      {
         /* we have a dual-knapsack constraint where we either need to choose one variable out of a subset (big
          * coefficients) of all or two variables of the rest
          *
          * e.g. 9x1 + 9x2 + 6x3 + 4x4 <= 19   <=>    9~x1 + 9~x2 + 6~x3 + 4~x4 >= 9
          *                                    <=>    2~x1 + 2~x2 +  ~x3 +  ~x4 >= 2
          *                                    <=>    2x1  +  2x2 +   x3 +   x4 <= 4
          *
          *      3x1 + 3x2 + 2x3 + 2x4 + 2x5 + 2x6 + x7 <= 12   <=>   3~x1 + 3~x2 + 2~x3 + 2~x4 + 2~x5 + 2~x6 + ~x7 >= 3
          *                                                     <=>   2~x1 + 2~x2 +  ~x3 +  ~x4 +  ~x5 +  ~x6 + ~x7 >= 2
          *                                                     <=>   2 x1 + 2 x2 +   x3 +   x4 +   x5 +   x6 +  x7 <= 7
          *
          */
         if( v > 0 && weights[nvars - 2] > 1 )
         {
            int ncoefchg = 0;

            /* reduce all bigger weights */
            for( w = 0; w < v; ++w )
            {
               if( weights[w] > 2 )
               {
                  consdataChgWeight(consdata, w, 2LL);
                  ++ncoefchg;
               }
               else
               {
                  assert(weights[0] == 2);
                  assert(weights[v - 1] == 2);
                  break;
               }
            }

            /* reduce all smaller weights */
            for( w = v; w < nvars; ++w )
            {
               if( weights[w] > 1 )
               {
                  consdataChgWeight(consdata, w, 1LL);
                  ++ncoefchg;
               }
            }
            assert(ncoefchg > 0);

            (*nchgcoefs) += ncoefchg;

            /* correct the capacity */
            consdata->capacity = (-2 + v * 2 + nvars - v); /*lint !e647*/
            assert(consdata->capacity > 0);
            assert(weights[0] <= consdata->capacity);
            assert(consdata->weightsum > consdata->capacity);
            /* reset the reductionsum */
            reductionsum = 0;
         }
         else if( v == 0 )
         {
            assert(weights[nvars - 2] == 1);
         }
      }
      else
      {
         SCIP_Longint minweight = weights[nvars - 1];
         SCIP_Longint newweight = dualcapacity - minweight;
         SCIP_Longint restsumweights = 0;
         SCIP_Longint sumcoef;
         SCIP_Bool sumcoefcase = FALSE;
         int startv = v;
         int end;
         int k;

         assert(weights[nvars - 1] + weights[nvars - 2] <= capacity);

         /* reduce big weights of pairs that exceed the dualcapacity
          *
          * e.g. 9x1 + 9x2 + 6x3 + 4x4 + 4x5 + 4x6 <= 27   <=>    9~x1 + 9~x2 + 6~x3 + 4~x4 + 4~x5 + 4~x6 >= 9
          *                                                <=>    9~x1 + 9~x2 + 5~x3 + 4~x4 + 4~x5 + 4~x6 >= 9
          *                                                <=>    9x1  + 9x2  + 5x3  + 4x4  + 4x5  + 4x6  <= 27
          */
         while( weights[v] > newweight )
         {
            reductionsum += (weights[v] - newweight);
            consdataChgWeight(consdata, v, newweight);
            ++v;
            assert(v < nvars);
         }
         (*nchgcoefs) += (v - startv);

         /* skip equal weights */
         while( weights[v] == newweight )
            ++v;

         if( v > 0 )
         {
            for( w = v; w < nvars; ++w )
               restsumweights += weights[w];
         }
         else
            restsumweights = consdata->weightsum;

         if( restsumweights < dualcapacity )
         {
            /* we found redundant variables, which does not influence the feasibility of any integral solution, e.g.
             *
             * +61x1  + 61x2  + 61x3  + 61x4  + 61x5  + 61x6  + 35x7  + 10x8 <= 350  <=>
             * +61~x1 + 61~x2 + 61~x3 + 61~x4 + 61~x5 + 61~x6 + 35~x7 + 10~x8 >= 61
             */
            if( startv == v )
            {
               /* remove redundant variables */
               for( w = nvars - 1; w >= v; --w )
               {
                  SCIP_CALL( delCoefPos(scip, cons, v) );
                  ++(*nchgcoefs);
               }

#ifndef NDEBUG
               /* each coefficients should exceed the dualcapacity by itself */
               for( ; w >= 0; --w )
                  assert(weights[w] == dualcapacity);
#endif
               /* for performance reasons we do not update the capacity(, i.e. reduce it by reductionsum) and directly
                * upgrade this constraint
                */
               SCIP_CALL( upgradeCons(scip, cons, ndelconss, naddconss) );
               assert(SCIPconsIsDeleted(cons));

               return SCIP_OKAY;
            }

            /* special case where we have three different coefficient types
             *
             * e.g. 9x1 + 9x2 + 6x3 + 6x4 + 4x5 + 4x6 <= 29   <=>    9~x1 + 9~x2 + 6~x3 + 6~x4 + 4~x5 + 4~x6 >= 9
             *                                                <=>    9~x1 + 9~x2 + 5~x3 + 5~x4 + 4~x5 + 4~x6 >= 9
             *                                                <=>    3~x1 + 3~x2 + 2~x3 + 2~x4 +  ~x5 +  ~x6 >= 3
             *                                                <=>    3x1  + 3x2  + 2x3  + 2x4  +   x5 +   x6 <= 9
             */
            if( weights[v] > 1 || (weights[startv] > (SCIP_Longint)nvars - v) || (startv > 0 && weights[0] == (SCIP_Longint)nvars - v + 1) )
            {
               SCIP_Longint newcap;

               /* adjust smallest coefficients, which all together do not exceed the dualcapacity */
               for( w = nvars - 1; w >= v; --w )
               {
                  if( weights[w] > 1 )
                  {
                     consdataChgWeight(consdata, w, 1LL);
                     ++(*nchgcoefs);
                  }
               }

               /* adjust middle sized coefficients, which when choosing also one small coefficients exceed the
                * dualcapacity
                */
               newweight = (SCIP_Longint)nvars - v;
               assert(newweight > 1);
               for( ; w >= startv; --w )
               {
                  if( weights[w] > newweight )
                  {
                     consdataChgWeight(consdata, w, newweight);
                     ++(*nchgcoefs);
                  }
                  else
                     assert(weights[w] == newweight);
               }

               /* adjust big sized coefficients, where each of them exceeds the dualcapacity by itself */
               ++newweight;
               assert(newweight > 2);
               for( ; w >= 0; --w )
               {
                  if( weights[w] > newweight )
                  {
                     consdataChgWeight(consdata, w, newweight);
                     ++(*nchgcoefs);
                  }
                  else
                     assert(weights[w] == newweight);
               }

               /* update the capacity */
               newcap = ((SCIP_Longint)startv - 1) * newweight + ((SCIP_Longint)v - startv) * (newweight - 1)  + ((SCIP_Longint)nvars - v);
               if( consdata->capacity > newcap )
               {
                  consdata->capacity = newcap;
                  ++(*nchgsides);
               }
               else
                  assert(consdata->capacity == newcap);
            }
            assert(weights[v] == 1 && (weights[startv] == (SCIP_Longint)nvars - v) && (startv == 0 || weights[0] == (SCIP_Longint)nvars - v + 1));

            /* the new dualcapacity should still be equal to the (nvars - v + 1) */
            assert(consdata->weightsum - consdata->capacity == (SCIP_Longint)nvars - v + 1);

            /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
             * weight must not be sorted by their index
             */
#ifndef NDEBUG
            for( w = nvars - 1; w > 0; --w )
               assert(weights[w] <= weights[w - 1]);
#endif
            return SCIP_OKAY;
         }

         /* check if all rear items have the same weight as the last one, so we cannot tighten the constraint further */
         end = nvars - 2;
         while( end >= 0 && weights[end] == weights[end + 1] )
         {
            assert(end >= v);
            --end;
         }

         if( v >= end )
            goto TERMINATE;

         end = nvars - 2;

         /* can we stop early, another special reduction case might exist */
         if( 2 * weights[end] > dualcapacity )
         {
            restsumweights = 0;

            /* determine capacity of the small items */
            for( w = end + 1; w < nvars; ++w )
               restsumweights += weights[w];

            if( restsumweights * 2 <= dualcapacity )
            {
               /* check for further posssible reductions in the middle */
               while( v < end && restsumweights + weights[v] >= dualcapacity )
                  ++v;

               if( v >= end )
                  goto TERMINATE;

               /* dualcapacity is even, we can set the middle weights to dualcapacity/2 */
               if( (dualcapacity & 1) == 0 )
               {
                  newweight = dualcapacity / 2;

                  /* set all middle coefficients */
                  for( ; v <= end; ++v )
                  {
                     if( weights[v] > newweight )
                     {
                        reductionsum += (weights[v] - newweight);
                        consdataChgWeight(consdata, v, newweight);
                        ++(*nchgcoefs);
                     }
                  }
               }
               /* dualcapacity is odd, we can set the middle weights to dualcapacity but therefor need to multiply all
                * other coefficients by 2
                */
               else
               {
                  /* correct the reductionsum */
                  reductionsum *= 2;

                  /* multiply big coefficients by 2 */
                  for( w = 0; w < v; ++w )
                  {
                     consdataChgWeight(consdata, w, weights[w] * 2);
                  }

                  newweight = dualcapacity;
                  /* set all middle coefficients */
                  for( ; v <= end; ++v )
                  {
                     reductionsum += (2 * weights[v] - newweight);
                     consdataChgWeight(consdata, v, newweight);
                  }

                  /* multiply small coefficients by 2 */
                  for( w = end + 1; w < nvars; ++w )
                  {
                     consdataChgWeight(consdata, w, weights[w] * 2);
                  }
                  (*nchgcoefs) += nvars;

                  dualcapacity *= 2;
                  consdata->capacity *= 2;
                  ++(*nchgsides);
               }
            }

            goto TERMINATE;
         }

         /* further reductions using the next possible coefficient sum
          *
          * e.g. 9x1 + 8x2 + 7x3 + 3x4 + x5 <= 19   <=>    9~x1 + 8~x2 + 7~x3 + 3~x4 + ~x5 >= 9
          *                                         <=>    9~x1 + 8~x2 + 6~x3 + 3~x4 + ~x5 >= 9
          *                                         <=>    9x1  + 8x2  + 6x3  + 3x4  + x5  <= 18
          */
         /* @todo loop for "k" can be extended, same coefficient when determine next sumcoef can be left out */
         for( k = 0; k < 4; ++k )
         {
            /* determine next minimal coefficient sum */
            switch( k )
            {
            case 0:
               sumcoef = weights[nvars - 1] + weights[nvars - 2];
               break;
            case 1:
               assert(nvars >= 3);
               sumcoef = weights[nvars - 1] + weights[nvars - 3];
               break;
            case 2:
               assert(nvars >= 4);
               if( weights[nvars - 1] + weights[nvars - 4] < weights[nvars - 2] + weights[nvars - 3] )
               {
                  sumcoefcase = TRUE;
                  sumcoef = weights[nvars - 1] + weights[nvars - 4];
               }
               else
               {
                  sumcoefcase = FALSE;
                  sumcoef = weights[nvars - 2] + weights[nvars - 3];
               }
               break;
            case 3:
               assert(nvars >= 5);
               if( sumcoefcase )
               {
                  sumcoef = MIN(weights[nvars - 1] + weights[nvars - 5], weights[nvars - 2] + weights[nvars - 3]);
               }
               else
               {
                  sumcoef = MIN(weights[nvars - 1] + weights[nvars - 4], weights[nvars - 1] + weights[nvars - 2] + weights[nvars - 3]);
               }
               break;
            default:
               return SCIP_ERROR;
            }

            /* tighten next coefficients that, pair with the current small coefficient, exceed the dualcapacity */
            minweight = weights[end];
            while( minweight <= sumcoef )
            {
               newweight = dualcapacity - minweight;
               startv = v;
               assert(v < nvars);

               /* @todo check for further reductions, when two times the minweight exceeds the dualcapacity */
               /* shrink big coefficients */
               while( weights[v] + minweight > dualcapacity && 2 * minweight <= dualcapacity )
               {
                  reductionsum += (weights[v] - newweight);
                  consdataChgWeight(consdata, v, newweight);
                  ++v;
                  assert(v < nvars);
               }
               (*nchgcoefs) += (v - startv);

               /* skip unchangable weights */
               while( weights[v] + minweight == dualcapacity )
               {
                  assert(v < nvars);
                  ++v;
               }

               --end;
               /* skip same end weights */
               while( end >= 0 && weights[end] == weights[end + 1] )
                  --end;

               if( v >= end )
                  goto TERMINATE;

               minweight = weights[end];
            }

            if( v >= end )
               goto TERMINATE;

            /* now check if a combination of small coefficients allows us to tighten big coefficients further */
            if( sumcoef < minweight )
            {
               minweight = sumcoef;
               newweight = dualcapacity - minweight;
               startv = v;
               assert(v < nvars);

               /* shrink big coefficients */
               while( weights[v] + minweight > dualcapacity && 2 * minweight <= dualcapacity )
               {
                  reductionsum += (weights[v] - newweight);
                  consdataChgWeight(consdata, v, newweight);
                  ++v;
                  assert(v < nvars);
               }
               (*nchgcoefs) += (v - startv);

               /* skip unchangable weights */
               while( weights[v] + minweight == dualcapacity )
               {
                  assert(v < nvars);
                  ++v;
               }
            }

            if( v >= end )
               goto TERMINATE;

            /* can we stop early, another special reduction case might exist */
            if( 2 * weights[end] > dualcapacity )
            {
               restsumweights = 0;

               /* determine capacity of the small items */
               for( w = end + 1; w < nvars; ++w )
                  restsumweights += weights[w];

               if( restsumweights * 2 <= dualcapacity )
               {
                  /* check for further posssible reductions in the middle */
                  while( v < end && restsumweights + weights[v] >= dualcapacity )
                     ++v;

                  if( v >= end )
                     goto TERMINATE;

                  /* dualcapacity is even, we can set the middle weights to dualcapacity/2 */
                  if( (dualcapacity & 1) == 0 )
                  {
                     newweight = dualcapacity / 2;

                     /* set all middle coefficients */
                     for( ; v <= end; ++v )
                     {
                        if( weights[v] > newweight )
                        {
                           reductionsum += (weights[v] - newweight);
                           consdataChgWeight(consdata, v, newweight);
                           ++(*nchgcoefs);
                        }
                     }
                  }
                  /* dualcapacity is odd, we can set the middle weights to dualcapacity but therefor need to multiply all
                   * other coefficients by 2
                   */
                  else
                  {
                     /* correct the reductionsum */
                     reductionsum *= 2;

                     /* multiply big coefficients by 2 */
                     for( w = 0; w < v; ++w )
                     {
                        consdataChgWeight(consdata, w, weights[w] * 2);
                     }

                     newweight = dualcapacity;
                     /* set all middle coefficients */
                     for( ; v <= end; ++v )
                     {
                        reductionsum += (2 * weights[v] - newweight);
                        consdataChgWeight(consdata, v, newweight);
                     }

                     /* multiply small coefficients by 2 */
                     for( w = end + 1; w < nvars; ++w )
                     {
                        consdataChgWeight(consdata, w, weights[w] * 2);
                     }
                     (*nchgcoefs) += nvars;

                     dualcapacity *= 2;
                     consdata->capacity *= 2;
                     ++(*nchgsides);
                  }
               }

               goto TERMINATE;
            }

            /* cannot tighten any further */
            if( 2 * sumcoef > dualcapacity )
               goto TERMINATE;
         }
      }
   }


 TERMINATE:
   /* correct capacity */
   if( reductionsum > 0 )
   {
      assert(v > 0);

      consdata->capacity -= reductionsum;
      ++(*nchgsides);

      assert(consdata->weightsum - dualcapacity == consdata->capacity);
   }
   assert(weights[0] <= consdata->capacity);

   /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal
    * weight must not be sorted by their index
    */
#ifndef NDEBUG
   for( w = nvars - 1; w > 0; --w )
      assert(weights[w] <= weights[w - 1]);
#endif

   if( oldnchgcoefs < *nchgcoefs )
   {
      assert(!SCIPconsIsDeleted(cons));

      /* it might be that we can divide the weights by their greatest common divisor */
      normalizeWeights(cons, nchgcoefs, nchgsides);
   }
   else
   {
      assert(oldnchgcoefs == *nchgcoefs);
      assert(oldnchgsides == *nchgsides);
   }

   return SCIP_OKAY;
}


/** fixes variables with weights bigger than the capacity and delete redundant constraints, also sort weights */
static
SCIP_RETCODE prepareCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nfixedvars,         /**< pointer to store the amount of fixed variables */
   int*                  ndelconss,          /**< pointer to store the amount of deleted constraints */
   int*                  nchgcoefs           /**< pointer to store the amount of changed coefficients */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSDATA* consdata;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* no variables left, then delete constraint */
   if( nvars == 0 )
   {
      assert(consdata->capacity >= 0);

      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);

      return SCIP_OKAY;
   }

   /* sort items */
   sortItems(consdata);

   vars = consdata->vars;
   weights = consdata->weights;
   capacity = consdata->capacity;
   v = 0;

   /* check for weights bigger than the capacity */
   while( v < nvars && weights[v] > capacity )
   {
      SCIP_CALL( SCIPfixVar(scip, vars[v], 0.0, &infeasible, &fixed) );
      assert(!infeasible);

      if( fixed )
         ++(*nfixedvars);

      ++v;
   }

   /* if we fixed at least one variable we need to delete them from the constraint */
   if( v > 0 )
   {
      if( v == nvars )
      {
         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*ndelconss);

         return SCIP_OKAY;
      }

      /* delete all position from back to front */
      for( --v; v >= 0; --v )
      {
         SCIP_CALL( delCoefPos(scip, cons, v) );
         ++(*nchgcoefs);
      }

      /* sort items again because of deletion */
      sortItems(consdata);
      assert(vars == consdata->vars);
      assert(weights == consdata->weights);
   }
   assert(consdata->sorted);
   assert(weights[0] <= capacity);

   if( !SCIPisHugeValue(scip, (SCIP_Real) capacity) && consdata->weightsum <= capacity )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
   }

   return SCIP_OKAY;
}


/** tries to simplify weights and delete redundant variables in knapsack a^Tx <= capacity
 *
 *  1. use the duality between a^Tx <= capacity   <=>   -a^T~x <= capacity - weightsum to tighten weights, e.g.
 *
 *     11x1 + 10x2 + 7x3 + 5x4 + 5x5 <= 25    <=>   -10~x1 - 10~x2 - 7~x3 - 5~x4 - 5~x5 <= -13
 *
 *     the above constraint can be changed to
 *
 *     -8~x1 - 8~x2 - 7~x3 - 5~x4 - 5~x5 <= -12   <=>   8x1 + 8x2 + 7x3 + 5x4 + 5x5 <= 20
 *
 *  2. if variables in a constraint do not affect the (in-)feasibility of the constraint, we can delete them, e.g.
 *
 *     7x1 + 6x2 + 5x3 + 5x4 + x5 + x6 <= 20 => x5 and x6 are redundant and can be removed
 *
 *  3. Tries to use gcd information an all but one weight to change this not-included weight and normalize the
 *     constraint further, e.g.
 *
 *     9x1 + 6x2 + 6x3 + 5x4 <= 13   =>   9x1 + 6x2 + 6x3 + 6x4 <= 12   =>   3x1 + 2x2 + 2x3 + 2x4 <= 4   =>   4x1 + 2x2 + 2x3 + 2x4 <= 4
 *                                                                                                        =>   2x1 + x2 + x3 + x4 <= 2
 *     9x1 + 6x2 + 6x3 + 7x4 <= 13   =>   9x1 + 6x2 + 6x3 + 6x4 <= 12   =>   see above
 */
static
SCIP_RETCODE simplifyInequalities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nfixedvars,         /**< pointer to store the amount of fixed variables */
   int*                  ndelconss,          /**< pointer to store the amount of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the amount of changed sides */
   int*                  naddconss,          /**< pointer to count number of added constraints */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSDATA* consdata;
   SCIP_Longint* weights;
   SCIP_Longint restweight;
   SCIP_Longint newweight;
   SCIP_Longint weight;
   SCIP_Longint oldgcd;
   SCIP_Longint rest;
   SCIP_Longint gcd;
   int oldnchgcoefs;
   int oldnchgsides;
   int candpos;
   int candpos2;
   int offsetv;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(naddconss != NULL);
   assert(cutoff != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *cutoff = FALSE;

   /* remove double enties and also combinations of active and negated variables */
   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   assert(consdata->merged);
   if( *cutoff )
      return SCIP_OKAY;

   assert(consdata->capacity >= 0);

   /* fix variables with big coefficients and remove redundant constraints, sort weights */
   SCIP_CALL( prepareCons(scip, cons, nfixedvars, ndelconss, nchgcoefs) );

   if( SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;

   if( !SCIPisHugeValue(scip, (SCIP_Real) consdata->capacity) )
   {
      /* 1. dual weights tightening */
      SCIP_CALL( dualWeightsTightening(scip, cons, ndelconss, nchgcoefs, nchgsides, naddconss) );

      if( SCIPconsIsDeleted(cons) )
         return SCIP_OKAY;
      /* 2. delete redundant variables */
      SCIP_CALL( detectRedundantVars(scip, cons, ndelconss, nchgcoefs, nchgsides, naddconss) );

      if( SCIPconsIsDeleted(cons) )
         return SCIP_OKAY;
   }

   weights = consdata->weights;
   nvars = consdata->nvars;

#ifndef NDEBUG
   /* constraint might not be sorted, but the weights are already sorted */
   for( v = nvars - 1; v > 0; --v )
      assert(weights[v] <= weights[v-1]);
#endif

   /* determine greatest common divisor */
   gcd = weights[nvars - 1];
   for( v = nvars - 2; v >= 0 && gcd > 1; --v )
   {
      gcd = SCIPcalcGreComDiv(gcd, weights[v]);
   }

   /* divide the constraint by their greatest common divisor */
   if( gcd >= 2 )
   {
      for( v = nvars - 1; v >= 0; --v )
      {
         consdataChgWeight(consdata, v, weights[v]/gcd);
      }
      (*nchgcoefs) += nvars;

      consdata->capacity /= gcd;
      (*nchgsides)++;
   }
   assert(consdata->nvars == nvars);

   /* weight should still be sorted, because the reduction preserves this, but corresponding variables with equal weight
    * must not be sorted by their index
    */
#ifndef NDEBUG
   for( v = nvars - 1; v > 0; --v )
      assert(weights[v] <= weights[v-1]);
#endif

   /* 3. start gcd procedure for all variables */
   do
   {
      SCIPdebug( oldnchgcoefs = *nchgcoefs; )
      SCIPdebug( oldnchgsides = *nchgsides; )

      vars = consdata->vars;
      weights = consdata->weights;
      nvars = consdata->nvars;

      /* stop if we have two coefficients which are one in absolute value */
      if( weights[nvars - 1] == 1 && weights[nvars - 2] == 1 )
         return SCIP_OKAY;

      v = 0;
      /* determine coefficients as big as the capacity, these we do not need to take into account when calculating the
       * gcd
       */
      while( weights[v] == consdata->capacity )
      {
         ++v;
         assert(v < nvars);
      }

      /* all but one variable are as big as the capacity, this is handled elsewhere */
      if( v == nvars - 1 )
         return SCIP_OKAY;

      offsetv = v;

      gcd = -1;
      candpos = -1;
      candpos2 = -1;

      /* calculate greatest common divisor over all integer and binary variables and determine the candidate where we might
       * change the coefficient
       */
      for( v = nvars - 1; v >= offsetv; --v )
      {
         weight = weights[v];
         assert(weight >= 1);

         oldgcd = gcd;

         if( gcd == -1 )
         {
            gcd = weights[v];
            assert(gcd >= 1);
         }
         else
         {
            /* calculate greatest common divisor for all variables */
            gcd = SCIPcalcGreComDiv(gcd, weight);
         }

         /* if the greatest commmon divisor has become 1, we might have found the possible coefficient to change or we
          * can terminate
          */
         if( gcd == 1 )
         {
            /* found candidate */
            if( candpos == -1 )
            {
               gcd = oldgcd;
               candpos = v;

               /* if both first coefficients have a gcd of 1, both are candidates for the coefficient change */
               if( v == nvars - 2 )
                  candpos2 = v + 1;
            }
            /* two different variables lead to a gcd of one, so we cannot change a coefficient */
            else
            {
               if( candpos == v + 1 && candpos2 == v + 2 )
               {
                  assert(candpos2 == nvars - 1);

                  /* take new candidates */
                  candpos = candpos2;

                  /* recalculate gcd from scratch */
                  gcd = weights[v+1];
                  assert(gcd >= 1);

                  /* calculate greatest common divisor for variables */
                  gcd = SCIPcalcGreComDiv(gcd, weights[v]);
                  if( gcd == 1 )
                     return SCIP_OKAY;
               }
               else
                  /* cannot determine a possible coefficient for reduction */
                  return SCIP_OKAY;
            }
         }
      }
      assert(gcd >= 2);

      /* we should have found one coefficient, that led to a gcd of 1, otherwise we could normalize the constraint
       * further
       */
      assert(((candpos >= offsetv) || (candpos == -1 && offsetv > 0)) && candpos < nvars);

      /* determine the remainder of the capacity and the gcd */
      rest = consdata->capacity % gcd;
      assert(rest >= 0);
      assert(rest < gcd);

      if( candpos == -1 )
      {
         /* we assume that the constraint was normalized */
         assert(rest > 0);

         /* replace old with new capacity */
         consdata->capacity -= rest;
         ++(*nchgsides);

         /* replace old big coefficients with new capacity */
         for( v = 0; v < offsetv; ++v )
         {
            consdataChgWeight(consdata, v, consdata->capacity);
         }

         *nchgcoefs += offsetv;
         goto CONTINUE;
      }

      /* determine the remainder of the coefficient candidate and the gcd */
      restweight = weights[candpos] % gcd;
      assert(restweight >= 1);
      assert(restweight < gcd);

      /* calculate new coefficient */
      if( restweight > rest )
         newweight = weights[candpos] - restweight + gcd;
      else
         newweight = weights[candpos] - restweight;

      assert(newweight == 0 || SCIPcalcGreComDiv(gcd, newweight) == gcd);

      SCIPdebugMsg(scip, "gcd = %" SCIP_LONGINT_FORMAT ", rest = %" SCIP_LONGINT_FORMAT ", restweight = %" SCIP_LONGINT_FORMAT "; possible new weight of variable <%s> %" SCIP_LONGINT_FORMAT ", possible new capacity %" SCIP_LONGINT_FORMAT ", offset of coefficients as big as capacity %d\n", gcd, rest, restweight, SCIPvarGetName(vars[candpos]), newweight, consdata->capacity - rest, offsetv);

      /* must not change weights and capacity if one variable would be removed and we have a big coefficient,
       * e.g., 11x1 + 6x2 + 6x3 + 5x4 <= 11 => gcd = 6, offsetv = 1 => newweight = 0, but we would lose x1 = 1 => x4 = 0
       */
      if( newweight == 0 && offsetv > 0 )
         return SCIP_OKAY;

      if( rest > 0 )
      {
         /* replace old with new capacity */
         consdata->capacity -= rest;
         ++(*nchgsides);

         /* replace old big coefficients with new capacity */
         for( v = 0; v < offsetv; ++v )
         {
            consdataChgWeight(consdata, v, consdata->capacity);
         }

         *nchgcoefs += offsetv;
      }

      if( newweight == 0 )
      {
         /* delete redundant coefficient */
         SCIP_CALL( delCoefPos(scip, cons, candpos) );
         assert(consdata->nvars == nvars - 1);
         --nvars;
      }
      else
      {
         /* replace old with new coefficient */
         consdataChgWeight(consdata, candpos, newweight);
      }
      ++(*nchgcoefs);

      assert(consdata->vars == vars);
      assert(consdata->nvars == nvars);
      assert(consdata->weights == weights);

   CONTINUE:
      /* now constraint can be normalized, dividing it by the gcd */
      for( v = nvars - 1; v >= 0; --v )
      {
         consdataChgWeight(consdata, v, weights[v]/gcd);
      }
      (*nchgcoefs) += nvars;

      consdata->capacity /= gcd;
      ++(*nchgsides);

      SCIPdebugPrintCons(scip, cons, NULL);

      SCIPdebugMsg(scip, "we did %d coefficient changes and %d side changes on constraint %s when applying one round of the gcd algorithm\n", *nchgcoefs - oldnchgcoefs, *nchgsides - oldnchgsides, SCIPconsGetName(cons));
   }
   while( nvars >= 2 );

   return SCIP_OKAY;
}


/** inserts an element into the list of binary zero implications */
static
SCIP_RETCODE insertZerolist(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 liftcands,          /**< array of the lifting candidates */
   int*                  nliftcands,         /**< number of lifting candidates */
   int**                 firstidxs,          /**< array of first zeroitems indices */
   SCIP_Longint**        zeroweightsums,     /**< array of sums of weights of the implied-to-zero items */
   int**                 zeroitems,          /**< pointer to zero items array */
   int**                 nextidxs,           /**< pointer to array of next zeroitems indeces */
   int*                  zeroitemssize,      /**< pointer to size of zero items array */
   int*                  nzeroitems,         /**< pointer to length of zero items array */
   int                   probindex,          /**< problem index of variable y in implication y == v -> x == 0 */
   SCIP_Bool             value,              /**< value v of variable y in implication */
   int                   knapsackidx,        /**< index of variable x in knapsack */
   SCIP_Longint          knapsackweight,     /**< weight of variable x in knapsack */
   SCIP_Bool*            memlimitreached     /**< pointer to store whether the memory limit was reached */
   )
{
   int nzeros;

   assert(liftcands != NULL);
   assert(liftcands[value] != NULL);
   assert(nliftcands != NULL);
   assert(firstidxs != NULL);
   assert(firstidxs[value] != NULL);
   assert(zeroweightsums != NULL);
   assert(zeroweightsums[value] != NULL);
   assert(zeroitems != NULL);
   assert(nextidxs != NULL);
   assert(zeroitemssize != NULL);
   assert(nzeroitems != NULL);
   assert(*nzeroitems <= *zeroitemssize);
   assert(0 <= probindex && probindex < SCIPgetNVars(scip) - SCIPgetNContVars(scip));
   assert(memlimitreached != NULL);

   nzeros = *nzeroitems;

   /* allocate enough memory */
   if( nzeros == *zeroitemssize )
   {
      /* we explicitly construct the complete implication graph where the knapsack variables are involved;
       * this can be too huge - abort on memory limit
       */
      if( *zeroitemssize >= MAX_ZEROITEMS_SIZE )
      {
         SCIPdebugMsg(scip, "memory limit of %d bytes reached in knapsack preprocessing - abort collecting zero items\n",
            *zeroitemssize);
         *memlimitreached = TRUE;
         return SCIP_OKAY;
      }
      *zeroitemssize *= 2;
      *zeroitemssize = MIN(*zeroitemssize, MAX_ZEROITEMS_SIZE);
      SCIP_CALL( SCIPreallocBufferArray(scip, zeroitems, *zeroitemssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, nextidxs, *zeroitemssize) );
   }
   assert(nzeros < *zeroitemssize);

   if( *memlimitreached )
      *memlimitreached = FALSE;

   /* insert element */
   (*zeroitems)[nzeros] = knapsackidx;
   (*nextidxs)[nzeros] = firstidxs[value][probindex];
   if( firstidxs[value][probindex] == 0 )
   {
      liftcands[value][nliftcands[value]] = probindex;
      ++nliftcands[value];
   }
   firstidxs[value][probindex] = nzeros;
   ++(*nzeroitems);
   zeroweightsums[value][probindex] += knapsackweight;

   return SCIP_OKAY;
}

#define MAX_CLIQUELENGTH 50
/** applies rule (3) of the weight tightening procedure, which can lift other variables into the knapsack:
 *  (3) for a clique C let C(xi == v) := C \ {j: xi == v -> xj == 0}),
 *      let cliqueweightsum(xi == v) := sum(W(C(xi == v)))
 *      if cliqueweightsum(xi == v) < capacity:
 *      - fixing variable xi to v would make the knapsack constraint redundant
 *      - the weight of the variable or its negation (depending on v) can be increased as long as it has the same
 *        redundancy effect:
 *          wi'       := capacity - cliqueweightsum(xi == v)
 *      this rule can also be applied to binary variables not in the knapsack!
 */
static
SCIP_RETCODE tightenWeightsLift(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   int nbinvars;
   int* liftcands[2];          /* binary variables that have at least one entry in zeroitems */
   int* firstidxs[2];          /* first index in zeroitems for each binary variable/value pair, or zero for empty list */
   SCIP_Longint* zeroweightsums[2]; /* sums of weights of the implied-to-zero items */
   int* zeroitems;             /* item number in knapsack that is implied to zero */
   int* nextidxs;              /* next index in zeroitems for the same binary variable, or zero for end of list */
   int zeroitemssize;
   int nzeroitems;
   SCIP_Bool* zeroiteminserted[2];
   SCIP_Bool memlimitreached;
   int nliftcands[2];
   SCIP_Bool* cliqueused;
   SCIP_Bool* itemremoved;
   SCIP_Longint maxcliqueweightsum;
   SCIP_VAR** addvars;
   SCIP_Longint* addweights;
   SCIP_Longint addweightsum;
   int nvars;
   int cliquenum;
   int naddvars;
   int val;
   int i;

   int* tmpindices;
   SCIP_Bool* tmpboolindices;
   int* tmpindices2;
   SCIP_Bool* tmpboolindices2;
   int* tmpindices3;
   SCIP_Bool* tmpboolindices3;
   int tmp;
   int tmp2;
   int tmp3;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(nchgcoefs != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL); /* we are in presolve, so no LP row exists */
   assert(consdata->weightsum > consdata->capacity); /* otherwise, the constraint is redundant */
   assert(consdata->nvars > 0);
   assert(consdata->merged);

   nvars = consdata->nvars;

   /* check if the knapsack has too many items/cliques for applying this costly method */
   if( (!consdata->cliquepartitioned && nvars > MAX_USECLIQUES_SIZE) || consdata->ncliques > MAX_USECLIQUES_SIZE )
      return SCIP_OKAY;

   /* sort items, s.t. the heaviest one is in the first position */
   sortItems(consdata);

   if( !consdata->cliquepartitioned && nvars > MAX_USECLIQUES_SIZE )
      return SCIP_OKAY;

   /* we have to consider all integral variables since even integer and implicit integer variables can have binary bounds */
   nbinvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   assert(nbinvars > 0);
   binvars = SCIPgetVars(scip);

   /* get conshdlrdata to use cleared memory */
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* allocate temporary memory for the list of implied to zero variables */
   zeroitemssize = MIN(nbinvars, MAX_ZEROITEMS_SIZE); /* initial size of zeroitems buffer */
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcands[0], nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftcands[1], nbinvars) );

   assert(conshdlrdata->ints1size > 0);
   assert(conshdlrdata->ints2size > 0);
   assert(conshdlrdata->longints1size > 0);
   assert(conshdlrdata->longints2size > 0);

   /* next if conditions should normally not be true, because it means that presolving has created more binary variables
    * than binary + integer variables existed at the presolving initialization method, but for example if you would
    * transform all integers into their binary representation then it maybe happens
    */
   if( conshdlrdata->ints1size < nbinvars )
   {
      int oldsize = conshdlrdata->ints1size;

      conshdlrdata->ints1size = nbinvars;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->ints1, oldsize, conshdlrdata->ints1size) );
      BMSclearMemoryArray(&(conshdlrdata->ints1[oldsize]), conshdlrdata->ints1size - oldsize); /*lint !e866*/
   }
   if( conshdlrdata->ints2size < nbinvars )
   {
      int oldsize = conshdlrdata->ints2size;

      conshdlrdata->ints2size = nbinvars;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->ints2, oldsize, conshdlrdata->ints2size) );
      BMSclearMemoryArray(&(conshdlrdata->ints2[oldsize]), conshdlrdata->ints2size - oldsize); /*lint !e866*/
   }
   if( conshdlrdata->longints1size < nbinvars )
   {
      int oldsize = conshdlrdata->longints1size;

      conshdlrdata->longints1size = nbinvars;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->longints1, oldsize, conshdlrdata->longints1size) );
      BMSclearMemoryArray(&(conshdlrdata->longints1[oldsize]), conshdlrdata->longints1size - oldsize); /*lint !e866*/
   }
   if( conshdlrdata->longints2size < nbinvars )
   {
      int oldsize = conshdlrdata->longints2size;

      conshdlrdata->longints2size = nbinvars;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->longints2, oldsize, conshdlrdata->longints2size) );
      BMSclearMemoryArray(&(conshdlrdata->longints2[oldsize]), conshdlrdata->longints2size - oldsize); /*lint !e866*/
   }

   firstidxs[0] = conshdlrdata->ints1;
   firstidxs[1] = conshdlrdata->ints2;
   zeroweightsums[0] = conshdlrdata->longints1;
   zeroweightsums[1] = conshdlrdata->longints2;

   /* check for cleared arrays, all entries are zero */
#ifndef NDEBUG
   for( tmp = nbinvars - 1; tmp >= 0; --tmp )
   {
      assert(firstidxs[0][tmp] == 0);
      assert(firstidxs[1][tmp] == 0);
      assert(zeroweightsums[0][tmp] == 0);
      assert(zeroweightsums[1][tmp] == 0);
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &zeroitems, zeroitemssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nextidxs, zeroitemssize) );

   zeroitems[0] = -1; /* dummy element */
   nextidxs[0] = -1;
   nzeroitems = 1;
   nliftcands[0] = 0;
   nliftcands[1] = 0;

   assert(conshdlrdata->bools1size > 0);
   assert(conshdlrdata->bools2size > 0);

   /* next if conditions should normally not be true, because it means that presolving has created more binary variables
    * than binary + integer variables existed at the presolving initialization method, but for example if you would
    * transform all integers into their binary representation then it maybe happens
    */
   if( conshdlrdata->bools1size < nbinvars )
   {
      int oldsize = conshdlrdata->bools1size;

      conshdlrdata->bools1size = nbinvars;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->bools1, oldsize, conshdlrdata->bools1size) );
      BMSclearMemoryArray(&(conshdlrdata->bools1[oldsize]), conshdlrdata->bools1size - oldsize); /*lint !e866*/
   }
   if( conshdlrdata->bools2size < nbinvars )
   {
      int oldsize = conshdlrdata->bools2size;

      conshdlrdata->bools2size = nbinvars;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->bools2, oldsize, conshdlrdata->bools2size) );
      BMSclearMemoryArray(&(conshdlrdata->bools2[oldsize]), conshdlrdata->bools2size - oldsize); /*lint !e866*/
   }

   zeroiteminserted[0] = conshdlrdata->bools1;
   zeroiteminserted[1] = conshdlrdata->bools2;

   /* check for cleared arrays, all entries are zero */
#ifndef NDEBUG
   for( tmp = nbinvars - 1; tmp >= 0; --tmp )
   {
      assert(zeroiteminserted[0][tmp] == 0);
      assert(zeroiteminserted[1][tmp] == 0);
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpboolindices, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices2, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpboolindices2, 2 * nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpindices3, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpboolindices3, consdata->nvars) );
   tmp2 = 0;
   tmp3 = 0;

   memlimitreached = FALSE;
   for( i = 0; i < consdata->nvars && !memlimitreached; ++i )
   {
      SCIP_CLIQUE** cliques;
      SCIP_VAR* var;
      SCIP_Longint weight;
      SCIP_Bool value;
      int varprobindex;
      int ncliques;
      int j;

      tmp = 0;

      /* get corresponding active problem variable */
      var = consdata->vars[i];
      weight = consdata->weights[i];
      value = TRUE;
      SCIP_CALL( SCIPvarGetProbvarBinary(&var, &value) );
      varprobindex = SCIPvarGetProbindex(var);
      assert(0 <= varprobindex && varprobindex < nbinvars);

      /* update the zeroweightsum */
      zeroweightsums[!value][varprobindex] += weight; /*lint !e514*/
      tmpboolindices3[tmp3] = !value;
      tmpindices3[tmp3] = varprobindex;
      ++tmp3;

      /* initialize the arrays of inserted zero items */
      /* first add the implications (~x == 1 -> x == 0) */
      {
         SCIP_Bool implvalue;
         int probindex;

         probindex = SCIPvarGetProbindex(var);
         assert(0 <= probindex && probindex < nbinvars);

         implvalue = !value;

         /* insert the item into the list of the implied variable/value */
         assert( !zeroiteminserted[implvalue][probindex] );

         if( firstidxs[implvalue][probindex] == 0 )
         {
            tmpboolindices2[tmp2] = implvalue;
            tmpindices2[tmp2] = probindex;
            ++tmp2;
         }
         SCIP_CALL( insertZerolist(scip, liftcands, nliftcands, firstidxs, zeroweightsums,
               &zeroitems, &nextidxs, &zeroitemssize, &nzeroitems, probindex, implvalue, i, weight,
               &memlimitreached) );
         zeroiteminserted[implvalue][probindex] = TRUE;
         tmpboolindices[tmp] = implvalue;
         tmpindices[tmp] = probindex;
         ++tmp;
      }

      /* get the cliques where the knapsack item is member of with value 1 */
      ncliques = SCIPvarGetNCliques(var, value);
      cliques = SCIPvarGetCliques(var, value);
      for( j = 0; j < ncliques && !memlimitreached; ++j )
      {
         SCIP_VAR** cliquevars;
         SCIP_Bool* cliquevalues;
         int ncliquevars;
         int k;

         ncliquevars = SCIPcliqueGetNVars(cliques[j]);

         /* discard big cliques */
         if( ncliquevars > MAX_CLIQUELENGTH )
            continue;

         cliquevars = SCIPcliqueGetVars(cliques[j]);
         cliquevalues = SCIPcliqueGetValues(cliques[j]);

         for( k = ncliquevars - 1; k >= 0; --k )
         {
            SCIP_Bool implvalue;
            int probindex;

            if( var == cliquevars[k] )
               continue;

            probindex = SCIPvarGetProbindex(cliquevars[k]);
            if( probindex == -1 )
               continue;

            assert(0 <= probindex && probindex < nbinvars);
            implvalue = cliquevalues[k];

            /* insert the item into the list of the clique variable/value */
            if( !zeroiteminserted[implvalue][probindex] )
            {
               if( firstidxs[implvalue][probindex] == 0 )
               {
                  tmpboolindices2[tmp2] = implvalue;
                  tmpindices2[tmp2] = probindex;
                  ++tmp2;
               }

               SCIP_CALL( insertZerolist(scip, liftcands, nliftcands, firstidxs, zeroweightsums,
                     &zeroitems, &nextidxs, &zeroitemssize, &nzeroitems, probindex, implvalue, i, weight,
                     &memlimitreached) );
               zeroiteminserted[implvalue][probindex] = TRUE;
               tmpboolindices[tmp] = implvalue;
               tmpindices[tmp] = probindex;
               ++tmp;

               if( memlimitreached )
                  break;
            }
         }
      }
      /* clear zeroiteminserted */
      for( --tmp; tmp >= 0; --tmp)
         zeroiteminserted[tmpboolindices[tmp]][tmpindices[tmp]] = FALSE;
   }
   SCIPfreeBufferArray(scip, &tmpboolindices);

   /* calculate the clique partition and the maximal sum of weights using the clique information */
   assert(consdata->sorted);
   SCIP_CALL( calcCliquepartition(scip, conshdlrdata, consdata, TRUE, FALSE) );

   assert(conshdlrdata->bools3size > 0);

   /* next if condition should normally not be true, because it means that presolving has created more binary variables
    * in one constraint than binary + integer variables existed in the whole problem at the presolving initialization
    * method, but for example if you would transform all integers into their binary representation then it maybe happens
    */
   if( conshdlrdata->bools3size < consdata->nvars )
   {
      int oldsize = conshdlrdata->bools3size;

      conshdlrdata->bools3size = consdata->nvars;;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->bools3, oldsize, conshdlrdata->bools3size) );
      BMSclearMemoryArray(&(conshdlrdata->bools3[oldsize]), conshdlrdata->bools3size - oldsize); /*lint !e866*/
   }

   cliqueused = conshdlrdata->bools3;

   /* check for cleared array, all entries are zero */
#ifndef NDEBUG
   for( tmp = consdata->nvars - 1; tmp >= 0; --tmp )
      assert(cliqueused[tmp] == 0);
#endif

   maxcliqueweightsum = 0;
   tmp = 0;

   /* calculates maximal weight of cliques */
   for( i = 0; i < consdata->nvars; ++i )
   {
      cliquenum = consdata->cliquepartition[i];
      assert(0 <= cliquenum && cliquenum < consdata->nvars);

      if( !cliqueused[cliquenum] )
      {
         maxcliqueweightsum += consdata->weights[i];
         cliqueused[cliquenum] = TRUE;
         tmpindices[tmp] = cliquenum;
         ++tmp;
      }
   }
   /* clear cliqueused */
   for( --tmp; tmp >= 0; --tmp)
      cliqueused[tmp] = FALSE;

   assert(conshdlrdata->bools4size > 0);

   /* next if condition should normally not be true, because it means that presolving has created more binary variables
    * in one constraint than binary + integer variables existed in the whole problem at the presolving initialization
    * method, but for example if you would transform all integers into their binary representation then it maybe happens
    */
   if( conshdlrdata->bools4size < consdata->nvars )
   {
      int oldsize = conshdlrdata->bools4size;

      conshdlrdata->bools4size = consdata->nvars;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->bools4, oldsize, conshdlrdata->bools4size) );
      BMSclearMemoryArray(&conshdlrdata->bools4[oldsize], conshdlrdata->bools4size - oldsize); /*lint !e866*/
   }

   itemremoved = conshdlrdata->bools4;

   /* check for cleared array, all entries are zero */
#ifndef NDEBUG
   for( tmp = consdata->nvars - 1; tmp >= 0; --tmp )
      assert(itemremoved[tmp] == 0);
#endif

   /* for each binary variable xi and each fixing v, calculate the cliqueweightsum and update the weight of the
    * variable in the knapsack (this is sequence-dependent because the new or modified weights have to be
    * included in subsequent cliqueweightsum calculations)
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &addvars, 2*nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &addweights, 2*nbinvars) );
   naddvars = 0;
   addweightsum = 0;
   for( val = 0; val < 2 && addweightsum < consdata->capacity; ++val )
   {
      for( i = 0; i < nliftcands[val] && addweightsum < consdata->capacity; ++i )
      {
         SCIP_Longint cliqueweightsum;
         int probindex;
         int idx;
         int j;

         tmp = 0;

         probindex = liftcands[val][i];
         assert(0 <= probindex && probindex < nbinvars);

         /* ignore empty zero lists and variables that cannot be lifted anyways */
         if( firstidxs[val][probindex] == 0
            || maxcliqueweightsum - zeroweightsums[val][probindex] + addweightsum >= consdata->capacity )
            continue;

         /* mark the items that are implied to zero by setting the current variable to the current value */
         for( idx = firstidxs[val][probindex]; idx != 0; idx = nextidxs[idx] )
         {
            assert(0 < idx && idx < nzeroitems);
            assert(0 <= zeroitems[idx] && zeroitems[idx] < consdata->nvars);
            itemremoved[zeroitems[idx]] = TRUE;
         }

         /* calculate the residual cliqueweight sum */
         cliqueweightsum = addweightsum; /* the previously added items are single-element cliques */
         for( j = 0; j < consdata->nvars; ++j )
         {
            cliquenum = consdata->cliquepartition[j];
            assert(0 <= cliquenum && cliquenum < consdata->nvars);
            if( !itemremoved[j] )
            {
               if( !cliqueused[cliquenum] )
               {
                  cliqueweightsum += consdata->weights[j];
                  cliqueused[cliquenum] = TRUE;
                  tmpindices[tmp] = cliquenum;
                  ++tmp;
               }

               if( cliqueweightsum >= consdata->capacity )
                  break;
            }
         }

         /* check if the weight of the variable/value can be increased */
         if( cliqueweightsum < consdata->capacity )
         {
            SCIP_VAR* var;
            SCIP_Longint weight;

            /* insert the variable (with value TRUE) in the list of additional items */
            assert(naddvars < 2*nbinvars);
            var = binvars[probindex];
            if( val == FALSE )
            {
               SCIP_CALL( SCIPgetNegatedVar(scip, var, &var) );
            }
            weight = consdata->capacity - cliqueweightsum;
            addvars[naddvars] = var;
            addweights[naddvars] = weight;
            addweightsum += weight;
            naddvars++;

            SCIPdebugMsg(scip, "knapsack constraint <%s>: adding lifted item %" SCIP_LONGINT_FORMAT "<%s>\n",
               SCIPconsGetName(cons), weight, SCIPvarGetName(var));
         }

         /* clear itemremoved */
         for( idx = firstidxs[val][probindex]; idx != 0; idx = nextidxs[idx] )
         {
            assert(0 < idx && idx < nzeroitems);
            assert(0 <= zeroitems[idx] && zeroitems[idx] < consdata->nvars);
            itemremoved[zeroitems[idx]] = FALSE;
         }
         /* clear cliqueused */
         for( --tmp; tmp >= 0; --tmp)
            cliqueused[tmpindices[tmp]] = FALSE;
      }
   }
   SCIPfreeBufferArray(scip, &tmpindices);

   /* clear part of zeroweightsums */
   for( --tmp3; tmp3 >= 0; --tmp3)
      zeroweightsums[tmpboolindices3[tmp3]][tmpindices3[tmp3]] = 0;

   /* clear rest of zeroweightsums and firstidxs */
   for( --tmp2; tmp2 >= 0; --tmp2)
   {
      zeroweightsums[tmpboolindices2[tmp2]][tmpindices2[tmp2]] = 0;
      firstidxs[tmpboolindices2[tmp2]][tmpindices2[tmp2]] = 0;
   }

   SCIPfreeBufferArray(scip, &tmpindices2);
   SCIPfreeBufferArray(scip, &tmpindices3);
   SCIPfreeBufferArray(scip, &tmpboolindices2);
   SCIPfreeBufferArray(scip, &tmpboolindices3);

   /* add all additional item weights */
   for( i = 0; i < naddvars; ++i )
   {
      SCIP_CALL( addCoef(scip, cons, addvars[i], addweights[i]) );
   }
   *nchgcoefs += naddvars;

   if( naddvars > 0 )
   {
      /* if new items were added, multiple entries of the same variable are possible and we have to clean up the constraint */
      SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &addweights);
   SCIPfreeBufferArray(scip, &addvars);
   SCIPfreeBufferArray(scip, &nextidxs);
   SCIPfreeBufferArray(scip, &zeroitems);
   SCIPfreeBufferArray(scip, &liftcands[1]);
   SCIPfreeBufferArray(scip, &liftcands[0]);

   return SCIP_OKAY;
}

/** tightens item weights and capacity in presolving:
 *  given a knapsack sum(wi*xi) <= capacity
 *  (1) let weightsum := sum(wi)
 *      if weightsum - wi < capacity:
 *      - not using item i would make the knapsack constraint redundant
 *      - wi and capacity can be changed to have the same redundancy effect and the same results for
 *        fixing xi to zero or one, but with a reduced wi and tightened capacity to tighten the LP relaxation
 *      - change coefficients:
 *          wi'       := weightsum - capacity
 *          capacity' := capacity - (wi - wi')
 *  (2) increase weights from front to back(sortation is necessary) if there is no space left for another weight
 *      - determine the four(can be adjusted) minimal weightsums of the knapsack, i.e. in increasing order
 *        weights[nvars - 1], weights[nvars - 2], MIN(weights[nvars - 3], weights[nvars - 1] + weights[nvars - 2]),
 *        MIN(MAX(weights[nvars - 3], weights[nvars - 1] + weights[nvars - 2]), weights[nvars - 4]), note that there
 *        can be multiple times the same weight, this can be improved
 *      - check if summing up a minimal weightsum with a big weight exceeds the capacity, then we can increase the big
 *        weight, to capacity - lastmininmalweightsum, e.g. :
 *        19x1 + 15x2 + 10x3 + 5x4 + 5x5 <= 19
 *         ->  minimal weightsums: 5, 5, 10, 10
 *         ->  15 + 5 > 19 => increase 15 to 19 - 0 = 19
 *         ->  10 + 10 > 19 => increase 10 to 19 - 5 = 14, resulting in
 *        19x1 + 19x2 + 14x3 + 5x4 + 5x5  <= 19
 *  (3) let W(C) be the maximal weight of clique C,
 *          cliqueweightsum := sum(W(C))
 *      if cliqueweightsum - W(C) < capacity:
 *      - not using any item of C would make the knapsack constraint redundant
 *      - weights wi, i in C, and capacity can be changed to have the same redundancy effect and the same results for
 *        fixing xi, i in C, to zero or one, but with a reduced wi and tightened capacity to tighten the LP relaxation
 *      - change coefficients:
 *          delta     := capacity - (cliqueweightsum - W(C))
 *          wi'       := max(wi - delta, 0)
 *          capacity' := capacity - delta
 *      This rule has to add the used cliques in order to ensure they are enforced - otherwise, the reduction might
 *      introduce infeasible solutions.
 *  (4) for a clique C let C(xi == v) := C \ {j: xi == v -> xj == 0}),
 *      let cliqueweightsum(xi == v) := sum(W(C(xi == v)))
 *      if cliqueweightsum(xi == v) < capacity:
 *      - fixing variable xi to v would make the knapsack constraint redundant
 *      - the weight of the variable or its negation (depending on v) can be increased as long as it has the same
 *        redundancy effect:
 *          wi'       := capacity - cliqueweightsum(xi == v)
 *      This rule can also be applied to binary variables not in the knapsack!
 *  (5) if min{w} + wi > capacity:
 *      - using item i would force to fix other items to zero
 *      - wi can be increased to the capacity
 */
static
SCIP_RETCODE tightenWeights(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   SCIP_PRESOLTIMING     presoltiming,       /**< current presolving timing */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides,          /**< pointer to count number of side changes */
   int*                  naddconss,          /**< pointer to count number of added constraints */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Longint* weights;
   SCIP_Longint sumcoef;
   SCIP_Longint capacity;
   SCIP_Longint newweight;
   SCIP_Longint maxweight;
   SCIP_Longint minweight;
   SCIP_Bool sumcoefcase = FALSE;
   int startpos;
   int backpos;
   int nvars;
   int pos;
   int k;
   int i;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(!SCIPconsIsModifiable(cons));

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL); /* we are in presolve, so no LP row exists */
   assert(consdata->onesweightsum == 0); /* all fixed variables should have been removed */
   assert(consdata->weightsum > consdata->capacity); /* otherwise, the constraint is redundant */
   assert(consdata->nvars > 0);

   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* apply rule (1) */
   if( (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 )
   {
      do
      {
         assert(consdata->merged);

         /* sort items, s.t. the heaviest one is in the first position */
         sortItems(consdata);

         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_Longint weight;

            weight = consdata->weights[i];
            if( consdata->weightsum - weight < consdata->capacity )
            {
               newweight = consdata->weightsum - consdata->capacity;
               consdataChgWeight(consdata, i, newweight);
               consdata->capacity -= (weight - newweight);
               (*nchgcoefs)++;
               (*nchgsides)++;
               assert(!consdata->sorted);
               SCIPdebugMsg(scip, "knapsack constraint <%s>: changed weight of <%s> from %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT ", capacity from %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT "\n",
                  SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), weight, newweight,
                  consdata->capacity + (weight-newweight), consdata->capacity);
            }
            else
               break;
         }
      }
      while( !consdata->sorted && consdata->weightsum > consdata->capacity );
   }

   /* check for redundancy */
   if( consdata->weightsum <= consdata->capacity )
      return SCIP_OKAY;

   pos = 0;
   while( pos < consdata->nvars && consdata->weights[pos] == consdata->capacity )
      ++pos;

   sumcoef = 0;
   weights = consdata->weights;
   nvars = consdata->nvars;
   capacity = consdata->capacity;

   if( (presoltiming & (SCIP_PRESOLTIMING_FAST | SCIP_PRESOLTIMING_MEDIUM)) != 0 &&
      pos < nvars && weights[pos] + weights[pos + 1] > capacity )
   {
      /* further reductions using the next possible coefficient sum
       *
       * e.g. 19x1 + 15x2 + 10x3 + 5x4 + 5x5 <= 19  <=>  19x1 + 19x2 + 14x3 + 5x4 + 5x5  <= 19
       */
      /* @todo loop for "k" can be extended, same coefficient when determine next sumcoef can be left out */
      for( k = 0; k < 4; ++k )
      {
         newweight = capacity - sumcoef;

         /* determine next minimal coefficient sum */
         switch( k )
         {
         case 0:
            sumcoef = weights[nvars - 1];
            backpos = nvars - 1;
            break;
         case 1:
            sumcoef = weights[nvars - 2];
            backpos = nvars - 2;
            break;
         case 2:
            if( weights[nvars - 3] < weights[nvars - 1] + weights[nvars - 2] )
            {
               sumcoefcase = TRUE;
               sumcoef = weights[nvars - 3];
               backpos = nvars - 3;
            }
            else
            {
               sumcoefcase = FALSE;
               sumcoef = weights[nvars - 1] + weights[nvars - 2];
               backpos = nvars - 2;
            }
            break;
         default:
            assert(k == 3);
            if( sumcoefcase )
            {
               if( weights[nvars - 4] < weights[nvars - 1] + weights[nvars - 2] )
               {
                  sumcoef = weights[nvars - 4];
                  backpos = nvars - 4;
               }
               else
               {
                  sumcoef = weights[nvars - 1] + weights[nvars - 2];
                  backpos = nvars - 2;
               }
            }
            else
            {
               sumcoef = weights[nvars - 3];
               backpos = nvars - 3;
            }
            break;
         }

         if( backpos <= pos )
            break;

         /* tighten next coefficients that, paired with the current small coefficient, exceed the capacity */
         maxweight = weights[pos];
         startpos = pos;
         while( 2 * maxweight > capacity && maxweight + sumcoef > capacity )
         {
            assert(newweight > weights[pos]);

            SCIPdebugMsg(scip, "in constraint <%s> changing weight %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT "\n",
               SCIPconsGetName(cons), maxweight, newweight);

            consdataChgWeight(consdata, pos, newweight);

            ++pos;
            assert(pos < nvars);

            maxweight = weights[pos];

            if( backpos <= pos )
               break;
         }
         (*nchgcoefs) += (pos - startpos);

         /* skip unchangable weights */
         while( pos < nvars && weights[pos] + sumcoef == capacity )
            ++pos;

         /* check special case were there is only one weight left to tighten
          *
          * e.g.  95x1 + 59x2 + 37x3 + 36x4 <= 95 (37 > 36)
          *
          *   =>  95x1 + 59x2 + 59x3 + 36x4 <= 95
          *
          *       197x1 + 120x2 + 77x3 + 10x4 <= 207 (here we cannot tighten the coefficient further)
          */
         if( pos + 1 == backpos && weights[pos] > sumcoef &&
            ((k == 0) || (k == 1 && weights[nvars - 1] + sumcoef + weights[pos] > capacity)) )
         {
            newweight = capacity - sumcoef;
            assert(newweight > weights[pos]);

            SCIPdebugMsg(scip, "in constraint <%s> changing weight %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT "\n",
               SCIPconsGetName(cons), maxweight, newweight);

            consdataChgWeight(consdata, pos, newweight);

            break;
         }

         if( backpos <= pos )
            break;
      }
   }

   /* apply rule (2) (don't apply, if the knapsack has too many items for applying this costly method) */
   if( (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 )
   {
      if( conshdlrdata->disaggregation && consdata->nvars - pos <= MAX_USECLIQUES_SIZE && consdata->nvars >= 2 &&
         pos > 0 && (SCIP_Longint)consdata->nvars - pos <= consdata->capacity &&
         consdata->weights[pos - 1] == consdata->capacity && (pos == consdata->nvars || consdata->weights[pos] == 1) )
      {
         SCIP_VAR** clqvars;
         SCIP_CONS* cliquecons;
         char name[SCIP_MAXSTRLEN];
         int* clqpart;
         int nclqvars;
         int nclq;
         int len;
         int c;
         int w;

         assert(!SCIPconsIsDeleted(cons));

         if( pos == consdata->nvars )
         {
            SCIPdebugMsg(scip, "upgrading knapsack constraint <%s> to a set-packing constraint", SCIPconsGetName(cons));

            SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, SCIPconsGetName(cons), pos, consdata->vars,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                  SCIPconsIsStickingAtNode(cons)) );

            SCIP_CALL( SCIPaddCons(scip, cliquecons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
            ++(*naddconss);

            /* delete old constraint */
            SCIP_CALL( SCIPdelCons(scip, cons) );
            ++(*ndelconss);

            return SCIP_OKAY;
         }

         len = consdata->nvars - pos;

         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &clqpart, len) );

         /* calculate clique partition */
         SCIP_CALL( SCIPcalcCliquePartition(scip, &(consdata->vars[pos]), len, clqpart, &nclq) );
         assert(nclq <= len);

#ifndef NDEBUG
         /* clique numbers must be at least as high as the index */
         for( w = 0; w < nclq; ++w )
            assert(clqpart[w] <= w);
#endif

         SCIPdebugMsg(scip, "Disaggregating knapsack constraint <%s> due to clique information.\n", SCIPconsGetName(cons));

         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &clqvars, pos + len - nclq + 1) );

         /* copy corresponding variables with big coefficients */
         for( w = pos - 1; w >= 0; --w )
            clqvars[w] = consdata->vars[w];

         /* create for each clique a set-packing constraint */
         for( c = 0; c < nclq; ++c )
         {
            nclqvars = pos;

            for( w = c; w < len; ++w )
            {
               if( clqpart[w] == c )
               {
                  assert(nclqvars < pos + len - nclq + 1);
                  clqvars[nclqvars] = consdata->vars[w + pos];
                  ++nclqvars;
               }
            }

            assert(nclqvars > 1);

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clq_%" SCIP_LONGINT_FORMAT "_%d", SCIPconsGetName(cons), consdata->capacity, c);
            SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, nclqvars, clqvars,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                  SCIPconsIsStickingAtNode(cons)) );
            SCIPdebugMsg(scip, " -> adding clique constraint: ");
            SCIPdebugPrintCons(scip, cliquecons, NULL);
            SCIP_CALL( SCIPaddCons(scip, cliquecons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
            ++(*naddconss);
         }

         /* delete old constraint */
         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*ndelconss);

         SCIPfreeBufferArray(scip, &clqvars);
         SCIPfreeBufferArray(scip, &clqpart);

         return SCIP_OKAY;
      }
      else if( consdata->nvars <= MAX_USECLIQUES_SIZE || (consdata->cliquepartitioned && consdata->ncliques <= MAX_USECLIQUES_SIZE) )
      {
         SCIP_Longint* maxcliqueweights;
         SCIP_Longint* newweightvals;
         int* newweightidxs;
         SCIP_Longint cliqueweightsum;

         SCIP_CALL( SCIPallocBufferArray(scip, &maxcliqueweights, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &newweightvals, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &newweightidxs, consdata->nvars) );

         /* repeat as long as changes have been applied */
         do
         {
            int ncliques;
            int cliquenum;
            SCIP_Bool zeroweights;

            assert(consdata->merged);

            /* sort items, s.t. the heaviest one is in the first position */
            sortItems(consdata);

            /* calculate a clique partition */
            SCIP_CALL( calcCliquepartition(scip, conshdlrdata, consdata, TRUE, FALSE) );

            /* if there are only single element cliques, rule (2) is equivalent to rule (1) */
            if( consdata->cliquepartition[consdata->nvars - 1] == consdata->nvars - 1 )
               break;

            /* calculate the maximal weight of the cliques and store the clique type */
            cliqueweightsum = 0;
            ncliques = 0;

            for( i = 0; i < consdata->nvars; ++i )
            {
               SCIP_Longint weight;

               cliquenum = consdata->cliquepartition[i];
               assert(0 <= cliquenum && cliquenum <= ncliques);

               weight = consdata->weights[i];
               assert(weight > 0);

               if( cliquenum == ncliques )
               {
                  maxcliqueweights[ncliques] = weight;
                  cliqueweightsum += weight;
                  ++ncliques;
               }

               assert(maxcliqueweights[cliquenum] >= weight);
            }

            /* apply rule on every clique */
            zeroweights = FALSE;
            for( i = 0; i < ncliques; ++i )
            {
               SCIP_Longint delta;

               delta = consdata->capacity - (cliqueweightsum - maxcliqueweights[i]);
               if( delta > 0 )
               {
                  SCIP_Longint newcapacity;
#ifndef NDEBUG
                  SCIP_Longint newmincliqueweight;
#endif
                  SCIP_Longint newminweightsuminclique;
                  SCIP_Bool forceclique;
                  int nnewweights;
                  int j;

                  SCIPdebugMsg(scip, "knapsack constraint <%s>: weights of clique %d (maxweight: %" SCIP_LONGINT_FORMAT ") can be tightened: cliqueweightsum=%" SCIP_LONGINT_FORMAT ", capacity=%" SCIP_LONGINT_FORMAT " -> delta: %" SCIP_LONGINT_FORMAT "\n",
                     SCIPconsGetName(cons), i, maxcliqueweights[i], cliqueweightsum, consdata->capacity, delta);
                  newcapacity = consdata->capacity - delta;
                  forceclique = FALSE;
                  nnewweights = 0;
#ifndef NDEBUG
                  newmincliqueweight = newcapacity + 1;
                  for( j = 0; j < i; ++j )
                     assert(consdata->cliquepartition[j] < i); /* no element j < i can be in clique i */
#endif
                  for( j = i; j < consdata->nvars; ++j )
                  {
                     if( consdata->cliquepartition[j] == i )
                     {
                        newweight = consdata->weights[j] - delta;
                        newweight = MAX(newweight, 0);

                        /* cache the new weight */
                        assert(nnewweights < consdata->nvars);
                        newweightvals[nnewweights] = newweight;
                        newweightidxs[nnewweights] = j;
                        nnewweights++;

#ifndef NDEBUG
                        assert(newweight <= newmincliqueweight); /* items are sorted by non-increasing weight! */
                        newmincliqueweight = newweight;
#endif
                     }
                  }

                  /* check if our clique information results out of this knapsack constraint and if so check if we would loose the clique information */
                  if( nnewweights > 1 )
                  {
#ifndef NDEBUG
                     j = newweightidxs[nnewweights - 2];
                     assert(0 <= j && j < consdata->nvars);
                     assert(consdata->cliquepartition[j] == i);
                     j = newweightidxs[nnewweights - 1];
                     assert(0 <= j && j < consdata->nvars);
                     assert(consdata->cliquepartition[j] == i);
#endif

                     newminweightsuminclique = newweightvals[nnewweights - 2];
                     newminweightsuminclique += newweightvals[nnewweights - 1];

                     /* check if these new two minimal weights both fit into the knapsack;
                      * if this is true, we have to add a clique constraint in order to enforce the clique
                      * (otherwise, the knapsack might have been one of the reasons for the clique, and the weight
                      * reduction might be infeasible, i.e., allows additional solutions)
                      */
                     if( newminweightsuminclique <= newcapacity )
                        forceclique = TRUE;
                  }

                  /* check if we really want to apply the change */
                  if( conshdlrdata->disaggregation || !forceclique )
                  {
                     SCIPdebugMsg(scip, " -> change capacity from %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT " (forceclique:%u)\n",
                        consdata->capacity, newcapacity, forceclique);
                     consdata->capacity = newcapacity;
                     (*nchgsides)++;

                     for( k = 0; k < nnewweights; ++k )
                     {
                        j = newweightidxs[k];
                        assert(0 <= j && j < consdata->nvars);
                        assert(consdata->cliquepartition[j] == i);

                        /* apply the weight change */
                        SCIPdebugMsg(scip, " -> change weight of <%s> from %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT "\n",
                           SCIPvarGetName(consdata->vars[j]), consdata->weights[j], newweightvals[k]);
                        consdataChgWeight(consdata, j, newweightvals[k]);
                        (*nchgcoefs)++;
                        assert(!consdata->sorted);
                        zeroweights = zeroweights || (newweightvals[k] == 0);
                     }
                     /* if before the weight update at least one pair of weights did not fit into the knapsack and now fits,
                      * we have to make sure, the clique is enforced - the clique might have been constructed partially from
                      * this constraint, and by reducing the weights, this clique information is not contained anymore in the
                      * knapsack constraint
                      */
                     if( forceclique )
                     {
                        SCIP_CONS* cliquecons;
                        char name[SCIP_MAXSTRLEN];
                        SCIP_VAR** cliquevars;

                        SCIP_CALL( SCIPallocBufferArray(scip, &cliquevars, nnewweights) );
                        for( k = 0; k < nnewweights; ++k )
                           cliquevars[k] = consdata->vars[newweightidxs[k]];

                        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_clq_%" SCIP_LONGINT_FORMAT "_%d", SCIPconsGetName(cons), consdata->capacity, i);
                        SCIP_CALL( SCIPcreateConsSetpack(scip, &cliquecons, name, nnewweights, cliquevars,
                              SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                              SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                              SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                              SCIPconsIsStickingAtNode(cons)) );
                        SCIPdebugMsg(scip, " -> adding clique constraint: ");
                        SCIPdebugPrintCons(scip, cliquecons, NULL);
                        SCIP_CALL( SCIPaddCons(scip, cliquecons) );
                        SCIP_CALL( SCIPreleaseCons(scip, &cliquecons) );
                        SCIPfreeBufferArray(scip, &cliquevars);
                        (*naddconss)++;
                     }
                  }
               }
            }
            if( zeroweights )
            {
               SCIP_CALL( removeZeroWeights(scip, cons) );
            }
         }
         while( !consdata->sorted && consdata->weightsum > consdata->capacity );

         /* free temporary memory */
         SCIPfreeBufferArray(scip, &newweightidxs);
         SCIPfreeBufferArray(scip, &newweightvals);
         SCIPfreeBufferArray(scip, &maxcliqueweights);

         /* check for redundancy */
         if( consdata->weightsum <= consdata->capacity )
            return SCIP_OKAY;
      }
   }

   /* apply rule (3) */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      SCIP_CALL( tightenWeightsLift(scip, cons, nchgcoefs, cutoff) );
   }

   /* check for redundancy */
   if( consdata->weightsum <= consdata->capacity )
      return SCIP_OKAY;

   if( (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 )
   {
      /* apply rule (4) (all but smallest weight) */
      assert(consdata->merged);
      sortItems(consdata);
      minweight = consdata->weights[consdata->nvars-1];
      for( i = 0; i < consdata->nvars-1; ++i )
      {
         SCIP_Longint weight;

         weight = consdata->weights[i];
         assert(weight >= minweight);
         if( minweight + weight > consdata->capacity )
         {
            if( weight < consdata->capacity )
            {
               SCIPdebugMsg(scip, "knapsack constraint <%s>: changed weight of <%s> from %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT "\n",
                  SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), weight, consdata->capacity);
               assert(consdata->sorted);
               consdataChgWeight(consdata, i, consdata->capacity); /* this does not destroy the weight order! */
               assert(i == 0 || consdata->weights[i-1] >= consdata->weights[i]);
               consdata->sorted = TRUE;
               (*nchgcoefs)++;
            }
         }
         else
            break;
      }

      /* apply rule (5) (smallest weight) */
      if( consdata->nvars >= 2 )
      {
         SCIP_Longint weight;

         minweight = consdata->weights[consdata->nvars-2];
         weight = consdata->weights[consdata->nvars-1];
         assert(minweight >= weight);
         if( minweight + weight > consdata->capacity && weight < consdata->capacity )
         {
            SCIPdebugMsg(scip, "knapsack constraint <%s>: changed weight of <%s> from %" SCIP_LONGINT_FORMAT " to %" SCIP_LONGINT_FORMAT "\n",
               SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[consdata->nvars-1]), weight, consdata->capacity);
            assert(consdata->sorted);
            consdataChgWeight(consdata, consdata->nvars-1, consdata->capacity); /* this does not destroy the weight order! */
            assert(minweight >= consdata->weights[consdata->nvars-1]);
            consdata->sorted = TRUE;
            (*nchgcoefs)++;
         }
      }
   }

   return SCIP_OKAY;
}


#ifdef SCIP_DEBUG
static
void printClique(
   SCIP_VAR**            cliquevars,
   int                   ncliquevars
   )
{
   int b;
   SCIPdebugMessage("adding new Clique: ");
   for( b = 0; b < ncliquevars; ++b )
      SCIPdebugPrintf("%s ", SCIPvarGetName(cliquevars[b]));
   SCIPdebugPrintf("\n");
}
#endif

/** adds negated cliques of the knapsack constraint to the global clique table */
static
SCIP_RETCODE addNegatedCliques(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< knapsack constraint */
   SCIP_Bool*const       cutoff,             /**< pointer to store whether the node can be cut off */
   int*const             nbdchgs             /**< pointer to count the number of performed bound changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** poscliquevars;
   SCIP_VAR** cliquevars;
   SCIP_Longint* maxweights;
   SCIP_Longint* gainweights;
   int* gaincliquepartition;
   SCIP_Bool* cliqueused;
   SCIP_Longint minactduetonegcliques;
   SCIP_Longint freecapacity;
   SCIP_Longint lastweight;
   SCIP_Longint beforelastweight;
   int nposcliquevars;
   int ncliquevars;
   int nvars;
   int nnegcliques;
   int lastcliqueused;
   int thisnbdchgs;
   int v;
   int w;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nbdchgs != NULL);

   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* check whether the cliques have already been added */
   if( consdata->cliquesadded || nvars == 0 )
      return SCIP_OKAY;

   /* make sure, the items are merged */
   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* make sure, items are sorted by non-increasing weight */
   sortItems(consdata);

   assert(consdata->merged);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   /* calculate a clique partition */
   SCIP_CALL( calcCliquepartition(scip, conshdlrdata, consdata, FALSE, TRUE) );
   nnegcliques = consdata->nnegcliques;

   /* if we have no negated cliques, stop */
   if( nnegcliques == nvars )
      return SCIP_OKAY;

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &poscliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gainweights, nvars) );
   BMSclearMemoryArray(gainweights, nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &gaincliquepartition, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxweights, nnegcliques) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliqueused, nnegcliques) );
   BMSclearMemoryArray(cliqueused, nnegcliques);

   nnegcliques = 0;
   minactduetonegcliques = 0;

   /* determine maximal weights for all negated cliques and calculate minimal weightsum due to negated cliques */
   for( v = 0; v < nvars; ++v )
   {
      assert(0 <= consdata->negcliquepartition[v] && consdata->negcliquepartition[v] <= nnegcliques);
      assert(consdata->weights[v] > 0);

      if( consdata->negcliquepartition[v] == nnegcliques )
      {
         nnegcliques++;
         maxweights[consdata->negcliquepartition[v]] = consdata->weights[v];
      }
      else
         minactduetonegcliques += consdata->weights[v];
   }

   nposcliquevars = 0;

   /* add cliques, using negated cliques information */
   if( minactduetonegcliques > 0 )
   {
      /* free capacity is the rest of not used capacity if the smallest amount of weights due to negated cliques are used */
      freecapacity = consdata->capacity - minactduetonegcliques;

      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMsg(scip, "Try to add negated cliques in knapsack constraint handler for constraint %s; capacity = %" SCIP_LONGINT_FORMAT ", minactivity(due to neg. cliques) = %" SCIP_LONGINT_FORMAT ", freecapacity = %" SCIP_LONGINT_FORMAT ".\n",
         SCIPconsGetName(cons), consdata->capacity, minactduetonegcliques, freecapacity);

      /* calculate possible gain by switching chosen items in negated cliques */
      for( v = 0; v < nvars; ++v )
      {
         if( !cliqueused[consdata->negcliquepartition[v]] )
         {
            cliqueused[consdata->negcliquepartition[v]] = TRUE;
            for( w = v + 1; w < nvars; ++w )
            {
               /* if we would take the biggest weight instead of another what would we gain, take weight[v] instead of
                * weight[w] (which are both in a negated clique) */
               if( consdata->negcliquepartition[v] == consdata->negcliquepartition[w] 
                  && consdata->weights[v] > consdata->weights[w] )
               {
                  poscliquevars[nposcliquevars] = consdata->vars[w];
                  gainweights[nposcliquevars] = maxweights[consdata->negcliquepartition[v]] - consdata->weights[w];
                  gaincliquepartition[nposcliquevars] = consdata->negcliquepartition[v];
                  ++nposcliquevars;
               }
            }
         }
      }

      /* try to create negated cliques */
      if( nposcliquevars > 0 )
      {
         /* sort possible gain per substitution of the clique members */
         SCIPsortDownLongPtrInt(gainweights,(void**) poscliquevars, gaincliquepartition, nposcliquevars);

         for( v = 0; v < nposcliquevars; ++v )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, poscliquevars[v], &cliquevars[0]) );
            ncliquevars = 1;
            lastweight = gainweights[v];
            beforelastweight = -1;
            lastcliqueused = gaincliquepartition[v];
            /* clear cliqueused to get an unused array */
            BMSclearMemoryArray(cliqueused, nnegcliques);
            cliqueused[gaincliquepartition[v]] = TRUE;

            /* taking bigger weights make the knapsack redundant so we will create cliques, only take items which are not
             * in the same negated clique and by taking two of them would exceed the free capacity */
            for( w = v + 1; w < nposcliquevars && !cliqueused[gaincliquepartition[w]] && gainweights[w] + lastweight > freecapacity; ++w )
            {
               beforelastweight = lastweight;
               lastweight = gainweights[w];
               lastcliqueused = gaincliquepartition[w];
               cliqueused[gaincliquepartition[w]] = TRUE;
               SCIP_CALL( SCIPgetNegatedVar(scip, poscliquevars[w], &cliquevars[ncliquevars]) );
               ++ncliquevars;
            }

            if( ncliquevars > 1 )
            {
               SCIPdebug( printClique(cliquevars, ncliquevars) );
               assert(beforelastweight > 0);
               /* add the clique to the clique table */
               /* this really happens, e.g., on enigma.mps from the short test set */
               SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, FALSE, cutoff, &thisnbdchgs) );
               if( *cutoff )
                  goto TERMINATE;
               *nbdchgs += thisnbdchgs;

               /* reset last used clique to get slightly different cliques */
               cliqueused[lastcliqueused] = FALSE;

               /* try to replace the last item in the clique by a different item to obtain a slightly different clique */
               for( ++w; w < nposcliquevars && !cliqueused[gaincliquepartition[w]] && beforelastweight + gainweights[w] > freecapacity; ++w )
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, poscliquevars[w], &cliquevars[ncliquevars - 1]) );
                  SCIPdebug( printClique(cliquevars, ncliquevars) );
                  SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, FALSE, cutoff, &thisnbdchgs) );
                  if( *cutoff )
                     goto TERMINATE;
                  *nbdchgs += thisnbdchgs;
               }
            }
         }
      }   
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cliqueused);
   SCIPfreeBufferArray(scip, &gaincliquepartition);
   SCIPfreeBufferArray(scip, &maxweights);
   SCIPfreeBufferArray(scip, &gainweights);
   SCIPfreeBufferArray(scip, &cliquevars);
   SCIPfreeBufferArray(scip, &poscliquevars);

   return SCIP_OKAY;
}

/** greedy clique detection by considering weights and capacity
 *
 *  greedily detects cliques by first sorting the items by decreasing weights (optional) and then collecting greedily
 *  1) neighboring items which exceed the capacity together => one clique
 *  2) looping through the remaining items and finding the largest set of preceding items to build a clique => possibly many more cliques
 */
static
SCIP_RETCODE greedyCliqueAlgorithm(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**            items,              /**< array of variable items */
   SCIP_Longint*         weights,            /**< weights of the items */
   int                   nitems,             /**< the number of items */
   SCIP_Longint          capacity,           /**< maximum free capacity of the knapsack */
   SCIP_Bool             sorteditems,        /**< are the items sorted by their weights nonincreasing? */
   SCIP_Real             cliqueextractfactor,/**< lower clique size limit for greedy clique extraction algorithm (relative to largest clique) */
   SCIP_Bool*const       cutoff,             /**< pointer to store whether the node can be cut off */
   int*const             nbdchgs             /**< pointer to count the number of performed bound changes */
   )
{
   SCIP_Longint lastweight;
   int ncliquevars;
   int i;
   int thisnbdchgs;

   if( nitems <= 1 )
      return SCIP_OKAY;

   /* sort possible gain per substitution of the clique members */
   if( ! sorteditems )
      SCIPsortDownLongPtr(weights,(void**) items, nitems);

   ncliquevars = 1;
   lastweight = weights[0];

   /* taking these two weights together violates the knapsack => include into clique */
   for( i = 1; i < nitems && weights[i] + lastweight > capacity; ++i )
   {
      lastweight = weights[i];
      ++ncliquevars;
   }

   if( ncliquevars > 1 )
   {
      SCIP_Longint compareweight;
      SCIP_VAR** cliquevars;
      int compareweightidx;
      int minclqsize;
      int nnzadded;

      /* add the clique to the clique table */
      SCIPdebug( printClique(items, ncliquevars) );
      SCIP_CALL( SCIPaddClique(scip, items, NULL, ncliquevars, FALSE, cutoff, &thisnbdchgs) );

      if( *cutoff )
         return SCIP_OKAY;

      *nbdchgs += thisnbdchgs;
      nnzadded = ncliquevars;

      /* no more cliques to be found (don't know if this can actually happen, since the knapsack could be replaced by a set-packing constraint)*/
      if( ncliquevars == nitems )
         return SCIP_OKAY;

      /* copy items in order into buffer array and deduce more cliques */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &cliquevars, items, ncliquevars) );

      /* try to replace the last item in the clique by a different item to obtain a slightly different clique */
      /* loop over remaining, smaller items and compare each item backwards against larger weights, starting with the second smallest weight */
      compareweightidx = ncliquevars - 2;
      assert(i == nitems || weights[i] + weights[ncliquevars - 1] <= capacity);

      /* determine minimum clique size for the following loop */
      minclqsize = (int)(cliqueextractfactor * ncliquevars);
      minclqsize = MAX(minclqsize, 2);

      /* loop over the remaining variables and the larger items of the first clique until we
       * find another clique or reach the size limit */
      while( compareweightidx >= 0 && i < nitems && ! (*cutoff)
            && ncliquevars >= minclqsize  /* stop at a given minimum clique size */
            && nnzadded <= 2 * nitems     /* stop if enough nonzeros were added to the cliquetable */
            )
      {
         compareweight = weights[compareweightidx];
         assert(compareweight > 0);

         /* include this item together with all items that have a weight at least as large as the compare weight in a clique */
         if( compareweight + weights[i] > capacity )
         {
            assert(compareweightidx == ncliquevars -2);
            cliquevars[ncliquevars - 1] = items[i];
            SCIPdebug( printClique(cliquevars, ncliquevars) );
            SCIP_CALL( SCIPaddClique(scip, cliquevars, NULL, ncliquevars, FALSE, cutoff, &thisnbdchgs) );

            nnzadded += ncliquevars;

            /* stop when there is a cutoff */
            if( ! (*cutoff) )
               *nbdchgs += thisnbdchgs;

            /* go to next smaller item */
            ++i;
         }
         else
         {
            /* choose a preceding, larger weight to compare small items against. Clique size is reduced by 1 simultaneously */
            compareweightidx--;
            ncliquevars --;
         }

      }

      SCIPfreeBufferArray(scip, &cliquevars);
   }

   return SCIP_OKAY;
}

/** adds cliques of the knapsack constraint to the global clique table */
static
SCIP_RETCODE addCliques(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< knapsack constraint */
   SCIP_Real             cliqueextractfactor,/**< lower clique size limit for greedy clique extraction algorithm (relative to largest clique) */
   SCIP_Bool*const       cutoff,             /**< pointer to store whether the node can be cut off */
   int*const             nbdchgs             /**< pointer to count the number of performed bound changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   SCIP_Longint minactduetonegcliques;
   SCIP_Longint freecapacity;
   int nnegcliques;
   int cliquenum;
   SCIP_VAR** poscliquevars;
   SCIP_Longint* gainweights;
   int nposcliquevars;
   SCIP_Longint* secondmaxweights;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nbdchgs != NULL);

   *cutoff = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* check whether the cliques have already been added */
   if( consdata->cliquesadded || nvars == 0 )
      return SCIP_OKAY;

   /* make sure, the items are merged */
   SCIP_CALL( mergeMultiples(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* make sure, the items are sorted by non-increasing weight */
   sortItems(consdata);

   assert(consdata->merged);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   /* calculate a clique partition */
   SCIP_CALL( calcCliquepartition(scip, conshdlrdata, consdata, FALSE, TRUE) );
   nnegcliques = consdata->nnegcliques;
   assert(nnegcliques <= nvars);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &poscliquevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &gainweights, nvars) );
   BMSclearMemoryArray(gainweights, nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &secondmaxweights, nnegcliques) );
   BMSclearMemoryArray(secondmaxweights, nnegcliques);

   minactduetonegcliques = 0;

   /* calculate minimal activity due to negated cliques, and determine second maximal weight in each clique */
   if( nnegcliques < nvars )
   {
      nnegcliques = 0;

      for( i = 0; i < nvars; ++i )
      {
         SCIP_Longint weight;

         cliquenum = consdata->negcliquepartition[i];
         assert(0 <= cliquenum && cliquenum <= nnegcliques);

         weight = consdata->weights[i];
         assert(weight > 0);

         if( cliquenum == nnegcliques )
            nnegcliques++;
         else
         {
            minactduetonegcliques += weight;
            if( secondmaxweights[cliquenum] == 0 )
               secondmaxweights[cliquenum] = weight;
         }
      }
   }

   /* add cliques, using negated cliques information */
   if( minactduetonegcliques > 0 )
   {
      /* free capacity is the rest of not used capacity if the smallest amount of weights due to negated cliques are used */
      freecapacity = consdata->capacity - minactduetonegcliques;

      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMsg(scip, "Try to add cliques in knapsack constraint handler for constraint %s; capacity = %" SCIP_LONGINT_FORMAT ", minactivity(due to neg. cliques) = %" SCIP_LONGINT_FORMAT ", freecapacity = %" SCIP_LONGINT_FORMAT ".\n",
         SCIPconsGetName(cons), consdata->capacity, minactduetonegcliques, freecapacity);

      /* create negated cliques out of negated cliques, if we do not take the smallest weight of a cliques ... */
      SCIP_CALL( addNegatedCliques(scip, cons, cutoff, nbdchgs ) );

      if( *cutoff )
         goto TERMINATE;

      nposcliquevars = 0;

      for( i = nvars - 1; i >= 0; --i )
      {
         /* if we would take the biggest weight instead of the second biggest */
         cliquenum = consdata->negcliquepartition[i];
         if( consdata->weights[i] > secondmaxweights[cliquenum] )
         {
            poscliquevars[nposcliquevars] = consdata->vars[i];
            gainweights[nposcliquevars] = consdata->weights[i] - secondmaxweights[cliquenum];
            ++nposcliquevars;
         }
      }

      /* use the gain weights and free capacity to derive greedily cliques */
      if( nposcliquevars > 1 )
      {
         SCIP_CALL( greedyCliqueAlgorithm(scip, poscliquevars, gainweights, nposcliquevars, freecapacity, FALSE, cliqueextractfactor, cutoff, nbdchgs) );

         if( *cutoff )
            goto TERMINATE;
      }
   }

   /* build cliques by using the items with the maximal weights */
   SCIP_CALL( greedyCliqueAlgorithm(scip, consdata->vars, consdata->weights, nvars, consdata->capacity, TRUE, cliqueextractfactor, cutoff, nbdchgs) );

   TERMINATE:
   /* free temporary memory and mark the constraint */
   SCIPfreeBufferArray(scip, &secondmaxweights);
   SCIPfreeBufferArray(scip, &gainweights);
   SCIPfreeBufferArray(scip, &poscliquevars);
   consdata->cliquesadded = TRUE;

   return SCIP_OKAY;
}


/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyKnapsackcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables and the 
 * same coefficients 
 */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqKnapsackcons)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   int i;

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   assert(consdata1->sorted);
   assert(consdata2->sorted);
#ifndef NDEBUG
   scip = (SCIP*)userptr; 
   assert(scip != NULL);
#endif

   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   for( i = consdata1->nvars - 1; i >= 0; --i )
   {
      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 || 
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 0); 

      /* tests if weights are equal too */  
      if( consdata1->weights[i] != consdata2->weights[i] )
         return FALSE;
   } 

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValKnapsackcons)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   SCIP_CONSDATA* consdata;
   int minidx;
   int mididx;
   int maxidx;

   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

#ifndef NDEBUG
   scip = (SCIP*)userptr; 
   assert(scip != NULL);
#endif

   /* sorts the constraints */
   sortItems(consdata);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   assert(minidx >= 0 && mididx >= 0 && maxidx >= 0);

   /* hash value depends on vectors of variable indices */
   return SCIPhashTwo(SCIPcombineFourInt(consdata->nvars, minidx, mididx, maxidx),
                      consdata->weights[0]);
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint 
 *  accordingly; in contrast to preprocessConstraintPairs(), it uses a hash table 
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the problem is infeasible */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(scip != NULL);
   assert(blkmem != NULL);
   assert(conss != NULL);
   assert(ndelconss != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = nconss;
   hashtablesize = MAX(hashtablesize, HASHSIZE_KNAPSACKCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyKnapsackcons, hashKeyEqKnapsackcons, hashKeyValKnapsackcons, (void*) scip) );

   /* check all constraints in the given set for redundancy */
   for( c = nconss - 1; c >= 0; --c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      consdata0 = SCIPconsGetData(cons0);
      assert(consdata0 != NULL);
      if( consdata0->nvars == 0 )
      {
         if( consdata0->capacity < 0 )
         {
            *cutoff = TRUE;
            goto TERMINATE;
         }
         else
         {
            SCIP_CALL( SCIPdelCons(scip, cons0) );
            ++(*ndelconss);
            continue;
         }
      }

      /* get constraint from current hash table with same variables and same weights as cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));

      if( cons1 != NULL )
      {
         SCIP_CONS* consstay;
         SCIP_CONS* consdel;
         SCIP_CONSDATA* consdata1;

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));

         /* constraint found: create a new constraint with same coefficients and best left and right hand side; 
          * delete old constraints afterwards
          */
         consdata1 = SCIPconsGetData(cons1);

         assert(consdata1 != NULL);
         assert(consdata0->nvars > 0 && consdata0->nvars == consdata1->nvars);

         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);
         assert(consdata0->weights[0] == consdata1->weights[0]);

         SCIPdebugMsg(scip, "knapsack constraints <%s> and <%s> with equal coefficients\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));

         /* check which constraint has to stay; */
         if( consdata0->capacity < consdata1->capacity )
         {
            consstay = cons0;
            consdel = cons1;

            /* exchange consdel with consstay in hashtable */
            SCIP_CALL( SCIPhashtableRemove(hashtable, (void*) consdel) );
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) consstay) );
         }
         else
         {
            consstay = cons1; 
            consdel = cons0; 
         }

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( SCIPupdateConsFlags(scip, consstay, consdel) );

         /* delete consdel */
         SCIP_CALL( SCIPdelCons(scip, consdel) );
         ++(*ndelconss);

         assert(SCIPconsIsActive(consstay));
      }
      else
      {
         /* no such constraint in current hash table: insert cons0 into hash table */  
         SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
      }
   }

 TERMINATE:
   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}


/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   int c;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(ndelconss != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(cons0 != NULL);
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   assert(consdata0->merged);

   /* sort the constraint */
   sortItems(consdata0);

   /* check constraint against all prior constraints */
   for( c = (consdata0->presolvedtiming == SCIP_PRESOLTIMING_EXHAUSTIVE ? firstchange : 0); c < chkind; ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_Bool iscons0incons1contained;
      SCIP_Bool iscons1incons0contained;
      SCIP_Real quotient;
      int v;
      int v0;
      int v1;

      cons1 = conss[c];
      assert(cons1 != NULL);
      if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* if both constraints didn't change since last pair processing, we can ignore the pair */
      if( consdata0->presolvedtiming >= SCIP_PRESOLTIMING_EXHAUSTIVE && consdata1->presolvedtiming >= SCIP_PRESOLTIMING_EXHAUSTIVE ) /*lint !e574*/
         continue;

      assert(consdata1->nvars >= 1);
      assert(consdata1->merged);

      /* sort the constraint */
      sortItems(consdata1);

      quotient = ((SCIP_Real) consdata0->capacity) / ((SCIP_Real) consdata1->capacity);

      if( consdata0->nvars > consdata1->nvars )
      {
         iscons0incons1contained = FALSE;
         iscons1incons0contained = TRUE;
         v = consdata1->nvars - 1;
      }
      else if( consdata0->nvars < consdata1->nvars )
      {
         iscons0incons1contained = TRUE;
         iscons1incons0contained = FALSE;
         v = consdata0->nvars - 1;
      }
      else
      {
         iscons0incons1contained = TRUE;
         iscons1incons0contained = TRUE;
         v = consdata0->nvars - 1;
      }

      SCIPdebugMsg(scip, "preprocess knapsack constraint pair <%s> and <%s>\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));

      /* check consdata0 against consdata1:
       * 1. if all variables var_i of cons1 are in cons0 and for each of these variables
       *    (consdata0->weights[i] / quotient) >= consdata1->weights[i] cons1 is redundant
       * 2. if all variables var_i of cons0 are in cons1 and for each of these variables
       *    (consdata0->weights[i] / quotient) <= consdata1->weights[i] cons0 is redundant
       */
      v0 = consdata0->nvars - 1;
      v1 = consdata1->nvars - 1;

      while( v >= 0 )
      {
         assert(iscons0incons1contained || iscons1incons0contained);

         /* now there are more variables in cons1 left */
         if( v1 > v0 )
         {
            iscons1incons0contained = FALSE;
            if( !iscons0incons1contained )
               break;
         }
         /* now there are more variables in cons0 left */
         else if( v1 < v0 )
         {
            iscons0incons1contained = FALSE;
            if( !iscons1incons0contained )
               break;
         }

         assert(v == v0 || v == v1);
	 assert(v0 >= 0);
	 assert(v1 >= 0);

         /* both variables are the same */
         if( consdata0->vars[v0] == consdata1->vars[v1] )
         {
            /* if cons1 is possible contained in cons0 (consdata0->weights[v0] / quotient) must be greater equals consdata1->weights[v1] */
            if( iscons1incons0contained && SCIPisLT(scip, ((SCIP_Real) consdata0->weights[v0]) / quotient, (SCIP_Real) consdata1->weights[v1]) )
            {
               iscons1incons0contained = FALSE;
               if( !iscons0incons1contained )
                  break;
            }
            /* if cons0 is possible contained in cons1 (consdata0->weight[v0] / quotient) must be less equals consdata1->weight[v1] */
            else if( iscons0incons1contained && SCIPisGT(scip, ((SCIP_Real) consdata0->weights[v0]) / quotient, (SCIP_Real) consdata1->weights[v1]) )
            {
               iscons0incons1contained = FALSE;
               if( !iscons1incons0contained )
                  break;
            }
            --v0;
            --v1;
	    --v;
         }
         else
         {
            /* both constraints have a variables which is not part of the other constraint, so stop */
            if( iscons0incons1contained && iscons1incons0contained )
            {
               iscons0incons1contained = FALSE;
               iscons1incons0contained = FALSE;
               break;
            }
            assert(iscons0incons1contained ? (v1 >= v0) : iscons1incons0contained);
            assert(iscons1incons0contained ? (v1 <= v0) : iscons0incons1contained);
            /* continue to the next variable */
            if( iscons0incons1contained )
               --v1;
            else
               --v0;
         }
      }
      /* neither one constraint was contained in another or we checked all variables of one constraint against the
       * other
       */
      assert(!iscons1incons0contained || !iscons0incons1contained || v0 == -1 || v1 == -1);

      if( iscons1incons0contained )
      {
         SCIPdebugMsg(scip, "knapsack constraint <%s> is redundant\n", SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons1, NULL);

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );

         SCIP_CALL( SCIPdelCons(scip, cons1) );
         ++(*ndelconss);
      }
      else if( iscons0incons1contained )
      {
         SCIPdebugMsg(scip, "knapsack constraint <%s> is redundant\n", SCIPconsGetName(cons0));
         SCIPdebugPrintCons(scip, cons0, NULL);

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( SCIPupdateConsFlags(scip, cons1, cons0) );

         SCIP_CALL( SCIPdelCons(scip, cons0) );
         ++(*ndelconss);
         break;
      }
   }

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
   SCIP_Bool violated;
   SCIP_Bool cutoff = FALSE;
   int maxncuts;
   int ncuts = 0;
   int i;

   *result = SCIP_FEASIBLE;

   SCIPdebugMsg(scip, "knapsack enforcement of %d/%d constraints for %s solution\n", nusefulconss, nconss,
         sol == NULL ? "LP" : "relaxation");

   /* get maximal number of cuts per round */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   maxncuts = (SCIPgetDepth(scip) == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   /* search for violated useful knapsack constraints */
   for( i = 0; i < nusefulconss && ncuts < maxncuts && ! cutoff; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, FALSE, FALSE, &violated) );
      if( violated )
      {
         /* add knapsack constraint as LP row to the relaxation */
         SCIP_CALL( addRelaxation(scip, conss[i], &cutoff) );
         ncuts++;
      }
   }

   /* as long as no violations were found, search for violated obsolete knapsack constraints */
   for( i = nusefulconss; i < nconss && ncuts == 0 && ! cutoff; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, FALSE, FALSE, &violated) );
      if( violated )
      {
         /* add knapsack constraint as LP row to the relaxation */
         SCIP_CALL( addRelaxation(scip, conss[i], &cutoff) );
         ncuts++;
      }
   }

   /* adjust the result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/*
 * Linear constraint upgrading
 */

/** creates and captures a knapsack constraint out of a linear inequality */
static
SCIP_RETCODE createNormalizedKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with inequality coefficients */
   SCIP_Real             lhs,                /**< left hand side of inequality */
   SCIP_Real             rhs,                /**< right hand side of inequality */
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
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? 
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_VAR** transvars;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Longint weight;
   int mult;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(SCIPisInfinity(scip, -lhs) != SCIPisInfinity(scip, rhs));

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

   /* if the right hand side is non-infinite, we have to negate all variables with negative coefficient;
    * otherwise, we have to negate all variables with positive coefficient and multiply the row with -1
    */
   if( SCIPisInfinity(scip, rhs) )
   {
      mult = -1;
      capacity = (SCIP_Longint)SCIPfeasFloor(scip, -lhs);
   }
   else
   {
      mult = +1;
      capacity = (SCIP_Longint)SCIPfeasFloor(scip, rhs);
   }

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPisFeasIntegral(scip, vals[v]));
      weight = mult * (SCIP_Longint)SCIPfeasFloor(scip, vals[v]);
      if( weight > 0 )
      {
         transvars[v] = vars[v];
         weights[v] = weight;
      }
      else
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
         weights[v] = -weight;
         capacity -= weight;
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   SCIP_CALL( SCIPcreateConsKnapsack(scip, cons, name, nvars, transvars, weights, capacity,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &transvars);

   return SCIP_OKAY;
}

/** tries to upgrade a linear constraint into a knapsack constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdKnapsack)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to a knapsack constraint
    * - all variables must be binary
    * - all coefficients must be integral
    * - exactly one of the sides must be infinite
    */
   upgrade = (nposbin + nnegbin + nposimplbin + nnegimplbin == nvars)
      && (ncoeffspone + ncoeffsnone + ncoeffspint + ncoeffsnint == nvars)
      && (SCIPisInfinity(scip, -lhs) != SCIPisInfinity(scip, rhs));

   if( upgrade )
   {
      SCIPdebugMsg(scip, "upgrading constraint <%s> to knapsack constraint\n", SCIPconsGetName(cons));

      /* create the knapsack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( createNormalizedKnapsack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
/**! [SnippetConsCopyKnapsack] */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyKnapsack)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}
/**! [SnippetConsCopyKnapsack] */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
/**! [SnippetConsFreeKnapsack] */
static
SCIP_DECL_CONSFREE(consFreeKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}
/**! [SnippetConsFreeKnapsack] */


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nvars;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* all variables which are of integral type can be binary; this can be checked via the method SCIPvarIsBinary(var) */
   nvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->reals1, nvars) );
   conshdlrdata->reals1size = nvars;

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->reals1, conshdlrdata->reals1size);
   conshdlrdata->reals1size = 0;

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nconss == 0 || conss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* all variables which are of integral type can be binary; this can be checked via the method SCIPvarIsBinary(var) */
   nvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->ints1, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->ints2, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->longints1, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->longints2, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->bools1, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->bools2, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->bools3, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->bools4, nvars) );

   conshdlrdata->ints1size = nvars;
   conshdlrdata->ints2size = nvars;
   conshdlrdata->longints1size = nvars;
   conshdlrdata->longints2size = nvars;
   conshdlrdata->bools1size = nvars;
   conshdlrdata->bools2size = nvars;
   conshdlrdata->bools3size = nvars;
   conshdlrdata->bools4size = nvars;

#ifdef WITH_CARDINALITY_UPGRADE
   conshdlrdata->upgradedcard = FALSE;
#endif

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   for( c = 0; c < nconss; ++c )
   {
      if( !SCIPconsIsDeleted(conss[c]) )
      {
         /* since we are not allowed to detect infeasibility in the exitpre stage, we dont give an infeasible pointer */
         SCIP_CALL( applyFixings(scip, conss[c], NULL) );
      }
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->ints1, conshdlrdata->ints1size);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->ints2, conshdlrdata->ints2size);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->longints1, conshdlrdata->longints1size);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->longints2, conshdlrdata->longints2size);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bools1, conshdlrdata->bools1size);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bools2, conshdlrdata->bools2size);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bools3, conshdlrdata->bools3size);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bools4, conshdlrdata->bools4size);

   conshdlrdata->ints1size = 0;
   conshdlrdata->ints2size = 0;
   conshdlrdata->longints1size = 0;
   conshdlrdata->longints2size = 0;
   conshdlrdata->bools1size = 0;
   conshdlrdata->bools2size = 0;
   conshdlrdata->bools3size = 0;
   conshdlrdata->bools4size = 0;

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* free knapsack constraint */
   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
/**! [SnippetConsTransKnapsack]*/
static
SCIP_DECL_CONSTRANS(consTransKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->nvars, sourcedata->vars, sourcedata->weights, sourcedata->capacity) ); 

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch events for variables */
   SCIP_CALL( catchEvents(scip, *targetcons, targetdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}
/**! [SnippetConsTransKnapsack]*/

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpKnapsack)
{  /*lint --e{715}*/
   int i;

   *infeasible = FALSE;

   for( i = 0; i < nconss && !(*infeasible); i++ )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i], infeasible) );
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool sepacardinality;
   SCIP_Bool cutoff;

   SCIP_Real loclowerbound;
   SCIP_Real glblowerbound;
   SCIP_Real cutoffbound;
   SCIP_Real maxbound;

   int depth;
   int nrounds;
   int sepafreq;
   int sepacardfreq;
   int ncuts;
   int maxsepacuts;
   int i;

   *result = SCIP_DIDNOTRUN;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   SCIPdebugMsg(scip, "knapsack separation of %d/%d constraints, round %d (max %d/%d)\n",
      nusefulconss, nconss, nrounds, conshdlrdata->maxroundsroot, conshdlrdata->maxrounds);

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* check, if we should additionally separate knapsack cuts */
   sepafreq = SCIPconshdlrGetSepaFreq(conshdlr);
   sepacardfreq = sepafreq * conshdlrdata->sepacardfreq;
   sepacardinality = (conshdlrdata->sepacardfreq >= 0)
      && ((sepacardfreq == 0 && depth == 0) || (sepacardfreq >= 1 && (depth % sepacardfreq == 0)));

   /* check dual bound to see if we want to produce knapsack cuts at this node */
   loclowerbound = SCIPgetLocalLowerbound(scip);
   glblowerbound = SCIPgetLowerbound(scip);
   cutoffbound = SCIPgetCutoffbound(scip);
   maxbound = glblowerbound + conshdlrdata->maxcardbounddist * (cutoffbound - glblowerbound);
   sepacardinality = sepacardinality && SCIPisLE(scip, loclowerbound, maxbound);
   sepacardinality = sepacardinality && (SCIPgetNLPBranchCands(scip) > 0);

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;
   cutoff = FALSE;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss && ncuts < maxsepacuts && !SCIPisStopped(scip); i++ )
   {
      SCIP_CALL( separateCons(scip, conss[i], NULL, sepacardinality, conshdlrdata->usegubs, &cutoff, &ncuts) );
   }

   /* adjust return value */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool sepacardinality;
   SCIP_Bool cutoff;

   int depth;
   int nrounds;
   int sepafreq;
   int sepacardfreq;
   int ncuts;
   int maxsepacuts;
   int i;

   *result = SCIP_DIDNOTRUN;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   SCIPdebugMsg(scip, "knapsack separation of %d/%d constraints, round %d (max %d/%d)\n",
      nusefulconss, nconss, nrounds, conshdlrdata->maxroundsroot, conshdlrdata->maxrounds);

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* check, if we should additionally separate knapsack cuts */
   sepafreq = SCIPconshdlrGetSepaFreq(conshdlr);
   sepacardfreq = sepafreq * conshdlrdata->sepacardfreq;
   sepacardinality = (conshdlrdata->sepacardfreq >= 0)
      && ((sepacardfreq == 0 && depth == 0) || (sepacardfreq >= 1 && (depth % sepacardfreq == 0)));

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;
   cutoff = FALSE;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss && ncuts < maxsepacuts && !SCIPisStopped(scip); i++ )
   {
      SCIP_CALL( separateCons(scip, conss[i], sol, sepacardinality, conshdlrdata->usegubs, &cutoff, &ncuts) );
   }

   /* adjust return value */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpKnapsack)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxKnapsack)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsKnapsack)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, TRUE, FALSE, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;  
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckKnapsack)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;

   for( i = 0; i < nconss && (*result == SCIP_FEASIBLE || completely); i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, printreason, &violated) );
      if( violated )
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   SCIP_Bool inpresolve;
   int nfixedvars;
   int i;

   cutoff = FALSE;
   nfixedvars = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   inpresolve = (SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);
   assert(!inpresolve || SCIPinProbing(scip));

   /* process useful constraints */
   for( i = 0; i < nmarkedconss && !cutoff; i++ )
   {
      /* do not propagate constraints with multi-aggregated variables, which should only happen in probing mode,
       * otherwise the multi-aggregation should be resolved
       */
      if( inpresolve && SCIPconsGetData(conss[i])->existmultaggr )
         continue;
#ifndef NDEBUG
      else
         assert(!(SCIPconsGetData(conss[i])->existmultaggr));
#endif

      SCIP_CALL( propagateCons(scip, conss[i], &cutoff, &redundant, &nfixedvars, conshdlrdata->negatedclique) );

      /* unmark the constraint to be propagated */
      SCIP_CALL( SCIPunmarkConsPropagate(scip, conss[i]) );

   }

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolKnapsack)
{  /*lint --e{574,715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   SCIP_Bool success;
   int oldnfixedvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnaddconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int firstchange;
   int c;
   SCIP_Bool newchanges;

   /* remember old preprocessing counters */
   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnaddconss = *naddconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;
   firstchange = INT_MAX;

   newchanges = (nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || nnewchgbds > 0 || nnewupgdconss > 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss && !SCIPisStopped(scip); c++ )
   {
      int thisnfixedvars;
      int thisnchgbds;

      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* update data structures */
      /* todo if UBTIGHTENED events were caught, we could move this block after the continue */
      if( newchanges || *nfixedvars > oldnfixedvars || *nchgbds > oldnchgbds )
      {
         SCIP_CALL( applyFixings(scip, cons, &cutoff) );
         if( cutoff )
            break;
      }

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->presolvedtiming = 0;
      else if( consdata->presolvedtiming >= presoltiming )
         continue;

      SCIPdebugMsg(scip, "presolving knapsack constraint <%s>\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);
      consdata->presolvedtiming = presoltiming;

      thisnfixedvars = *nfixedvars;
      thisnchgbds = *nchgbds;

      /* merge constraint, so propagation works better */
      SCIP_CALL( mergeMultiples(scip, cons, &cutoff) );
      if( cutoff )
         break;

      /* add cliques in the knapsack to the clique table */
      if( (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 )
      {
         SCIP_CALL( addCliques(scip, cons, conshdlrdata->cliqueextractfactor, &cutoff, nchgbds) );
         if( cutoff )
            break;
      }

      /* propagate constraint */
      if( presoltiming < SCIP_PRESOLTIMING_EXHAUSTIVE )
      {
         SCIP_CALL( propagateCons(scip, cons, &cutoff, &redundant, nfixedvars, (presoltiming & SCIP_PRESOLTIMING_MEDIUM)) );

         if( cutoff )
            break;
         if( redundant )
         {
            (*ndelconss)++;
            continue;
         }
      }

      /* remove again all fixed variables, if further fixings were found */
      if( *nfixedvars > thisnfixedvars || *nchgbds > thisnchgbds )
      {
         SCIP_CALL( applyFixings(scip, cons, &cutoff) );
         if( cutoff )
            break;

         thisnfixedvars = *nfixedvars;
      }

      if( !SCIPconsIsModifiable(cons) )
      {
         /* check again for redundancy (applyFixings() might have decreased weightsum due to fixed-to-zero vars) */
         if( consdata->weightsum <= consdata->capacity )
         {
            SCIPdebugMsg(scip, " -> knapsack constraint <%s> is redundant: weightsum=%" SCIP_LONGINT_FORMAT ", capacity=%" SCIP_LONGINT_FORMAT "\n",
               SCIPconsGetName(cons), consdata->weightsum, consdata->capacity);
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            continue;
         }

         /* divide weights by their greatest common divisor */
         normalizeWeights(cons, nchgcoefs, nchgsides);

         /* try to simplify inequalities */
         if( conshdlrdata->simplifyinequalities && (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 )
         {
            SCIP_CALL( simplifyInequalities(scip, cons, nfixedvars, ndelconss, nchgcoefs, nchgsides, naddconss, &cutoff) );
            if( cutoff )
               break;

            if( SCIPconsIsDeleted(cons) )
               continue;

            /* remove again all fixed variables, if further fixings were found */
            if( *nfixedvars > thisnfixedvars )
            {
               SCIP_CALL(applyFixings(scip, cons, &cutoff));
               if( cutoff )
                  break;
            }
         }

         /* tighten capacity and weights */
         SCIP_CALL( tightenWeights(scip, cons, presoltiming, nchgcoefs, nchgsides, naddconss, ndelconss, &cutoff) );
         if( cutoff )
            break;

         if( SCIPconsIsActive(cons) )
         {
            if( conshdlrdata->dualpresolving && SCIPallowDualReds(scip) && (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 )
            {
               /* in case the knapsack constraints is independent of everything else, solve the knapsack and apply the
                * dual reduction
                */
               SCIP_CALL( dualPresolving(scip, cons, nchgbds, ndelconss, &redundant) );
               if( redundant )
                  continue;
            }

            /* check if knapsack constraint is parallel to objective function */
            SCIP_CALL( checkParallelObjective(scip, cons, conshdlrdata) );
         }
      }
      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->presolvedtiming != SCIP_PRESOLTIMING_EXHAUSTIVE )
         firstchange = c;
   }

   /* preprocess pairs of knapsack constraints */
   if( !cutoff && conshdlrdata->presolusehashing && (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 )
   {
      /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
      SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &cutoff, ndelconss) );
   }

   if( (*ndelconss != oldndelconss) || (*nchgsides != oldnchgsides) || (*nchgcoefs != oldnchgcoefs) || (*naddconss != oldnaddconss) )
      success = TRUE;
   else
      success = FALSE;

   if( !cutoff && firstchange < nconss && conshdlrdata->presolpairwise && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      SCIP_Longint npaircomparisons;

      npaircomparisons = 0;
      oldndelconss = *ndelconss;
      oldnchgsides = *nchgsides;
      oldnchgcoefs = *nchgcoefs;

      for( c = firstchange; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
      {
         cons = conss[c];
         if( !SCIPconsIsActive(cons) || SCIPconsIsModifiable(cons) )
            continue;

         npaircomparisons += ((SCIPconsGetData(cons)->presolvedtiming < SCIP_PRESOLTIMING_EXHAUSTIVE) ? (SCIP_Longint) c : ((SCIP_Longint) c - (SCIP_Longint) firstchange));

         SCIP_CALL( preprocessConstraintPairs(scip, conss, firstchange, c, ndelconss) );

         if( npaircomparisons > NMINCOMPARISONS )
         {
            if( (*ndelconss != oldndelconss) || (*nchgsides != oldnchgsides) || (*nchgcoefs != oldnchgcoefs) )
               success = TRUE;
            if( ((SCIP_Real) (*ndelconss - oldndelconss) + ((SCIP_Real) (*nchgsides - oldnchgsides))/2.0 +
                  ((SCIP_Real) (*nchgcoefs - oldnchgcoefs))/10.0) / ((SCIP_Real) npaircomparisons) < MINGAINPERNMINCOMPARISONS )
               break;
            oldndelconss = *ndelconss;
            oldnchgsides = *nchgsides;
            oldnchgcoefs = *nchgcoefs;
            npaircomparisons = 0;
         }
      }
   }
#ifdef WITH_CARDINALITY_UPGRADE
   /* @todo upgrade to cardinality constraints: the code below relies on disabling the checking of the knapsack
    * constraint in the original problem, because the upgrade ensures that at most the given number of continuous
    * variables has a nonzero value, but not that the binary variables corresponding to the continuous variables with
    * value zero are set to zero as well. This can cause problems if the user accesses the values of the binary
    * variables (as the MIPLIB solution checker does), or the transformed problem is freed and the original problem
    * (possibly with some user modifications) is re-optimized. Until there is a way to force the binary variables to 0
    * as well, we better keep this code disabled. */
   /* upgrade to cardinality constraints - only try to upgrade towards the end of presolving, since the process below is quite expensive */
   if ( ! cutoff && conshdlrdata->upgdcardinality && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 && SCIPisPresolveFinished(scip) && ! conshdlrdata->upgradedcard )
   {
      SCIP_HASHMAP* varhash;
      SCIP_VAR** cardvars;
      SCIP_Real* cardweights;
      int noldupgdconss;
      int nscipvars;
      int makeupgrade;

      noldupgdconss = *nupgdconss;
      nscipvars = SCIPgetNVars(scip);
      SCIP_CALL( SCIPallocClearBufferArray(scip, &cardvars, nscipvars) );
      SCIP_CALL( SCIPallocClearBufferArray(scip, &cardweights, nscipvars) );

      /* set up hash map */
      SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip), nscipvars) );

      /* We loop through all cardinality constraints twice:
       * - First, determine for each binary variable the number of cardinality constraints that can be upgraded to a
       *   knapsack constraint and contain this variable; this number has to coincide with the number of variable up
       *   locks; otherwise it would be infeasible to delete the knapsack constraints after the constraint update.
       * - Second, upgrade knapsack constraints to cardinality constraints. */
      for (makeupgrade = 0; makeupgrade < 2; ++makeupgrade)
      {
         for (c = nconss-1; c >= 0 && ! SCIPisStopped(scip); --c)
         {
            SCIP_CONS* cardcons;
            SCIP_VAR** vars;
            SCIP_Longint* weights;
            int nvars;
            int v;

            cons = conss[c];
            assert( cons != NULL );
            consdata = SCIPconsGetData(cons);
            assert( consdata != NULL );

            nvars = consdata->nvars;
            vars = consdata->vars;
            weights = consdata->weights;

            /* Check, whether linear knapsack can be upgraded to a cardinality constraint:
             * - all variables must be binary (always true)
             * - all coefficients must be 1.0
             * - the right hand side must be smaller than nvars
             */
            if ( consdata->capacity >= nvars )
               continue;

            /* the weights are sorted: check first and last weight */
            assert( consdata->sorted );
            if ( weights[0] != 1 || weights[nvars-1] != 1 )
               continue;

            /* check whether all variables are of the form 0 <= x_v <= u_v y_v for y_v \in \{0,1\} and zero objective */
            for (v = 0; v < nvars; ++v)
            {
               SCIP_BOUNDTYPE* impltypes;
               SCIP_Real* implbounds;
               SCIP_VAR** implvars;
               SCIP_VAR* var;
               int nimpls;
               int j;

               var = consdata->vars[v];
               assert( var != NULL );
               assert( SCIPvarIsBinary(var) );

               /* ignore non-active variables */
               if ( ! SCIPvarIsActive(var) )
                  break;

               /* be sure that implication variable has zero objective */
               if ( ! SCIPisZero(scip, SCIPvarGetObj(var)) )
                  break;

               nimpls = SCIPvarGetNImpls(var, FALSE);
               implvars = SCIPvarGetImplVars(var, FALSE);
               implbounds = SCIPvarGetImplBounds(var, FALSE);
               impltypes = SCIPvarGetImplTypes(var, FALSE);

               for (j = 0; j < nimpls; ++j)
               {
                  /* be sure that continuous variable is fixed to 0 */
                  if ( impltypes[j] != SCIP_BOUNDTYPE_UPPER )
                     continue;

                  /* cannot currently deal with nonzero fixings */
                  if ( ! SCIPisZero(scip, implbounds[j]) )
                     continue;

                  /* number of down locks should be one */
                  if ( SCIPvarGetNLocksDown(vars[v]) != 1 )
                     continue;

                  cardvars[v] = implvars[j];
                  cardweights[v] = (SCIP_Real) v;

                  break;
               }

               /* found no variable upper bound candidate -> exit */
               if ( j >= nimpls )
                  break;
            }

            /* did not find fitting variable upper bound for some variable -> exit */
            if ( v < nvars )
               break;

            /* save number of knapsack constraints that can be upgraded to a cardinality constraint,
             * in which the binary variable is involved in */
            if ( makeupgrade == 0 )
            {
               for (v = 0; v < nvars; ++v)
               {
                  if ( SCIPhashmapExists(varhash, vars[v]) )
                  {
                     int image;

                     image = (int) (size_t) SCIPhashmapGetImage(varhash, vars[v]);
                     SCIP_CALL( SCIPhashmapSetImage(varhash, vars[v], (void*) (size_t) (image + 1)) );/*lint !e776*/
                     assert( image + 1 == (int) (size_t) SCIPhashmapGetImage(varhash, vars[v]) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPhashmapInsert(varhash, vars[v], (void*) (size_t) 1) );/*lint !e571*/
                     assert( 1 == (int) (size_t) SCIPhashmapGetImage(varhash, vars[v]) );
                     assert( SCIPhashmapExists(varhash, vars[v]) );
                  }
               }
            }
            else
            {
               SCIP_CONS* origcons;

               /* for each variable: check whether the number of cardinality constraints that can be upgraded to a
                * knapsack constraint coincides with the number of variable up locks */
               for (v = 0; v < nvars; ++v)
               {
                  assert( SCIPhashmapExists(varhash, vars[v]) );
                  if ( SCIPvarGetNLocksUp(vars[v]) != (int) (size_t) SCIPhashmapGetImage(varhash, vars[v]) )
                     break;
               }
               if ( v < nvars )
                  break;

               /* store that we have upgraded */
               conshdlrdata->upgradedcard = TRUE;

               /* at this point we found suitable variable upper bounds */
               SCIPdebugMessage("Upgrading knapsack constraint <%s> to cardinality constraint ...\n", SCIPconsGetName(cons));

               /* create cardinality constraint */
               assert( ! SCIPconsIsModifiable(cons) );
               SCIP_CALL( SCIPcreateConsCardinality(scip, &cardcons, SCIPconsGetName(cons), nvars, cardvars, (int) consdata->capacity, vars, cardweights,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, cons, NULL);
               SCIPinfoMessage(scip, NULL, "\n");
               SCIPprintCons(scip, cardcons, NULL);
               SCIPinfoMessage(scip, NULL, "\n");
#endif
               SCIP_CALL( SCIPaddCons(scip, cardcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cardcons) );
               ++(*nupgdconss);

               /* delete oknapsack constraint */
               SCIP_CALL( SCIPdelCons(scip, cons) );
               ++(*ndelconss);

               /* We need to disable the original knapsack constraint, since it might happen that the binary variables
                * are 1 although the continuous variables are 0. Thus, the knapsack constraint might be violated,
                * although the cardinality constraint is satisfied. */
               origcons = SCIPfindOrigCons(scip, SCIPconsGetName(cons));
               assert( origcons != NULL );
               SCIP_CALL( SCIPsetConsChecked(scip, origcons, FALSE) );

               for (v = 0; v < nvars; ++v)
               {
                  int image;

                  assert ( SCIPhashmapExists(varhash, vars[v]) );
                  image = (int) (size_t) SCIPhashmapGetImage(varhash, vars[v]);
                  SCIP_CALL( SCIPhashmapSetImage(varhash, vars[v], (void*) (size_t) (image - 1)) );
                  assert( image - 1 == (int) (size_t) SCIPhashmapGetImage(varhash, vars[v]) );
               }
            }
         }
      }
      SCIPhashmapFree(&varhash);
      SCIPfreeBufferArray(scip, &cardweights);
      SCIPfreeBufferArray(scip, &cardvars);

      if ( *nupgdconss > noldupgdconss )
         success = TRUE;
   }
#endif

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( success || *nfixedvars > oldnfixedvars || *nchgbds > oldnchgbds )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Longint capsum;
   int i;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if we fixed a binary variable to one (due to negated clique) */
   if( inferinfo >= 0 && SCIPvarGetLbLocal(infervar) > 0.5 )
   {
      for( i = 0; i < consdata->nvars; ++i )
      {
         if( SCIPvarGetIndex(consdata->vars[i]) == inferinfo ) 
         {
            assert( SCIPgetVarUbAtIndex(scip, consdata->vars[i], bdchgidx, FALSE) < 0.5 );
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
            break;
         }
      }
      assert(i < consdata->nvars);
   }
   else
   {
      /* according to negated cliques the minweightsum and all variables which are fixed to one which led to a fixing of
       * another negated clique variable to one, the inferinfo was chosen to be the negative of the position in the 
       * knapsack constraint, see one above call of SCIPinferBinvarCons
       */
      if( inferinfo < 0 )
         capsum = 0;
      else
      {
         /* locate the inference variable and calculate the capacity that has to be used up to conclude infervar == 0;
          * inferinfo stores the position of the inference variable (but maybe the variables were resorted)
          */
         if( inferinfo < consdata->nvars && consdata->vars[inferinfo] == infervar )
            capsum = consdata->weights[inferinfo];
         else
         {
            for( i = 0; i < consdata->nvars && consdata->vars[i] != infervar; ++i )
            {}
            assert(i < consdata->nvars);
            capsum = consdata->weights[i];
         }
      }

      /* add fixed-to-one variables up to the point, that their weight plus the weight of the conflict variable exceeds
       * the capacity
       */
      if( capsum <= consdata->capacity )
      {
         for( i = 0; i < consdata->nvars; i++ )
         {
            if( SCIPgetVarLbAtIndex(scip, consdata->vars[i], bdchgidx, FALSE) > 0.5 )
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[i]) );
               capsum += consdata->weights[i];
               if( capsum > consdata->capacity )
                  break;
            }
         }
      }
   }

   /* NOTE: It might be the case that capsum < consdata->capacity. This is due the fact that the fixing of the variable
    *       to zero can included negated clique information. A negated clique means, that at most one of the clique
    *       variables can be zero. These information can be used to compute a minimum activity of the constraint and
    *       used to fix variables to zero.
    *
    *       Even if capsum < consdata->capacity we still reported a complete reason since the minimum activity is based
    *       on global variable bounds. It might even be the case that we reported to many variables which are fixed to
    *       one.
    */
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
/**! [SnippetConsLockKnapsack] */
static
SCIP_DECL_CONSLOCK(consLockKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
   }

   return SCIP_OKAY;
}
/**! [SnippetConsLockKnapsack] */


/** variable deletion method of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelvarsKnapsack)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   if( nconss > 0 )
   {
      SCIP_CALL( performVarDeletions(scip, conshdlr, conss, nconss) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; ++i )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, " ");
      SCIPinfoMessage(scip, file, "%+" SCIP_LONGINT_FORMAT, consdata->weights[i]);
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[i], TRUE) );
   }
   SCIPinfoMessage(scip, file, " <= %" SCIP_LONGINT_FORMAT "", consdata->capacity);

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyKnapsack)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_Longint* weights;
   SCIP_Real* coefs;
   const char* consname;
   int nvars;
   int v;

   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsKnapsack(sourcescip, sourcecons);
   nvars = SCIPgetNVarsKnapsack(sourcescip, sourcecons);
   weights = SCIPgetWeightsKnapsack(sourcescip, sourcecons);

   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   for( v = 0; v < nvars; ++v )
      coefs[v] = (SCIP_Real) weights[v];

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* copy the logic using the linear constraint copy method */
   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, nvars, sourcevars, coefs,
         -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(sourcescip, sourcecons), varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL);

   SCIPfreeBufferArray(scip, &coefs);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseKnapsack)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_Longint weight;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   char* endptr;
   int nread;
   int nvars;
   int varssize;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = TRUE;

   nvars = 0;
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,    varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, varssize) );

   while( *str != '\0' )
   {
      /* try to parse coefficient, and stop if not successful (probably reached <=) */
      if( sscanf(str, "%" SCIP_LONGINT_FORMAT "%n", &weight, &nread) < 1 )
         break;

      str += nread;

      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;

      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, str, &var, &endptr) );
      if( var == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
         *success = FALSE;
         break;
      }

      str = endptr;

      /* store weight and variable */
      if( varssize <= nvars )
      {
         varssize = SCIPcalcMemGrowSize(scip, varssize+1);
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars,    varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &weights, varssize) );
      }

      vars[nvars]    = var;
      weights[nvars] = weight;
      ++nvars;

      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;
   }

   if( *success )
   {
      if( strncmp(str, "<= ", 3) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '<= ' at begin of '%s'\n", str);
         *success = FALSE;
      }
      else
      {
         str += 3;
      }
   }

   if( *success )
   {
      if( sscanf(str, "%" SCIP_LONGINT_FORMAT, &capacity) != 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "error parsing capacity from '%s'\n", str);
         *success = FALSE;
      }
      else
      {
	 SCIP_CALL( SCIPcreateConsKnapsack(scip, cons, name, nvars, vars, weights, capacity,
	       initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      }
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &weights);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Event handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecKnapsack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventdata != NULL);
   assert(eventdata->cons != NULL);

   consdata = SCIPconsGetData(eventdata->cons);
   assert(consdata != NULL);

   switch( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      consdata->onesweightsum += eventdata->weight;
      consdata->presolvedtiming = 0;
      SCIP_CALL( SCIPmarkConsPropagate(scip, eventdata->cons) );
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      consdata->onesweightsum -= eventdata->weight;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      consdata->presolvedtiming = 0;
      SCIP_CALL( SCIPmarkConsPropagate(scip, eventdata->cons) );
      break;
   case SCIP_EVENTTYPE_VARFIXED:  /* the variable should be removed from the constraint in presolving */
      if( !consdata->existmultaggr )
      {
         SCIP_VAR* var;
         var = SCIPeventGetVar(event);
         assert(var != NULL);

         /* if the variable was aggregated or multiaggregated, we must signal to propagation that we are no longer merged */
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
         {
            consdata->existmultaggr = TRUE;
            consdata->merged = FALSE;
         }
         else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED ||
            (SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED && SCIPvarGetStatus(SCIPvarGetNegatedVar(var)) == SCIP_VARSTATUS_AGGREGATED) )
            consdata->merged = FALSE;

      }
      /*lint -fallthrough*/
   case SCIP_EVENTTYPE_IMPLADDED: /* further preprocessing might be possible due to additional implications */
      consdata->presolvedtiming = 0;
      break;
   case SCIP_EVENTTYPE_VARDELETED:
      consdata->varsdeleted = TRUE;
      break;
   default:
      SCIPerrorMessage("invalid event type %x\n", SCIPeventGetType(event));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for knapsack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create knapsack constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   conshdlrdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(conshdlrdata->eventhdlr), EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecKnapsack, eventhdlrdata) );

   /* get event handler for bound change events */
   if( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for knapsack constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpKnapsack, consEnfopsKnapsack, consCheckKnapsack, consLockKnapsack,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyKnapsack, consCopyKnapsack) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteKnapsack) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsKnapsack) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitKnapsack) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreKnapsack) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolKnapsack) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeKnapsack) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsKnapsack) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsKnapsack) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitKnapsack) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreKnapsack) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpKnapsack) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseKnapsack) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolKnapsack,CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintKnapsack) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropKnapsack, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropKnapsack) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpKnapsack, consSepasolKnapsack, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransKnapsack) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxKnapsack) );

   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint to knapsack constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }

   /* add knapsack constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/sepacardfreq",
         "multiplier on separation frequency, how often knapsack cuts are separated (-1: never, 0: only at root)",
         &conshdlrdata->sepacardfreq, TRUE, DEFAULT_SEPACARDFREQ, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/maxcardbounddist",
         "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for separating knapsack cuts",
         &conshdlrdata->maxcardbounddist, TRUE, DEFAULT_MAXCARDBOUNDDIST, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/cliqueextractfactor",
         "lower clique size limit for greedy clique extraction algorithm (relative to largest clique)",
         &conshdlrdata->cliqueextractfactor, TRUE, DEFAULT_CLIQUEEXTRACTFACTOR, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/disaggregation",
         "should disaggregation of knapsack constraints be allowed in preprocessing?",
         &conshdlrdata->disaggregation, TRUE, DEFAULT_DISAGGREGATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/simplifyinequalities",
         "should presolving try to simplify knapsacks",
         &conshdlrdata->simplifyinequalities, TRUE, DEFAULT_SIMPLIFYINEQUALITIES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/negatedclique",
         "should negated clique information be used in solving process",
         &conshdlrdata->negatedclique, TRUE, DEFAULT_NEGATEDCLIQUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance", 
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/dualpresolving",
         "should dual presolving steps be performed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/usegubs",
         "should GUB information be used for separation?",
         &conshdlrdata->usegubs, TRUE, DEFAULT_USEGUBS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/detectcutoffbound",
         "should presolving try to detect constraints parallel to the objective function defining an upper bound and prevent these constraints from entering the LP?",
         &conshdlrdata->detectcutoffbound, TRUE, DEFAULT_DETECTCUTOFFBOUND, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/detectlowerbound",
         "should presolving try to detect constraints parallel to the objective function defining a lower bound and prevent these constraints from entering the LP?",
         &conshdlrdata->detectlowerbound, TRUE, DEFAULT_DETECTLOWERBOUND, NULL, NULL) );
    SCIP_CALL( SCIPaddBoolParam(scip,
          "constraints/" CONSHDLR_NAME "/updatecliquepartitions",
          "should clique partition information be updated when old partition seems outdated?",
          &conshdlrdata->updatecliquepartitions, TRUE, DEFAULT_UPDATECLIQUEPARTITIONS, NULL, NULL) );
    SCIP_CALL( SCIPaddRealParam(scip,
          "constraints/" CONSHDLR_NAME "/clqpartupdatefac",
          "factor on the growth of global cliques to decide when to update a previous "
          "(negated) clique partition (used only if updatecliquepartitions is set to TRUE)",
          &conshdlrdata->clqpartupdatefac, TRUE, DEFAULT_CLQPARTUPDATEFAC, 1.0, 10.0, NULL, NULL) );
#ifdef WITH_CARDINALITY_UPGRADE
    SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/upgdcardinality",
         "if TRUE then try to update knapsack constraints to cardinality constraints",
         &conshdlrdata->upgdcardinality, TRUE, DEFAULT_UPGDCARDINALITY, NULL, NULL) );
#endif
   return SCIP_OKAY;
}

/** creates and captures a knapsack constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
/**! [SnippetConsCreationKnapsack] */
SCIP_RETCODE SCIPcreateConsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of items in the knapsack */
   SCIP_VAR**            vars,               /**< array with item variables */
   SCIP_Longint*         weights,            /**< array with item weights */
   SCIP_Longint          capacity,           /**< capacity of knapsack (right hand side of inequality) */
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
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the knapsack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("knapsack constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars, weights, capacity) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   /* catch events for variables */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( catchEvents(scip, *cons, consdata, conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}
/**! [SnippetConsCreationKnapsack] */

/** creates and captures a knapsack constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsKnapsack(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsKnapsack() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of items in the knapsack */
   SCIP_VAR**            vars,               /**< array with item variables */
   SCIP_Longint*         weights,            /**< array with item weights */
   SCIP_Longint          capacity            /**< capacity of knapsack */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsKnapsack(scip, cons, name, nvars, vars, weights, capacity,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds new item to knapsack constraint */
SCIP_RETCODE SCIPaddCoefKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< item variable */
   SCIP_Longint          weight              /**< item weight */
   )
{
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( addCoef(scip, cons, var, weight) );

   return SCIP_OKAY;
}

/** gets the capacity of the knapsack constraint */
SCIP_Longint SCIPgetCapacityKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
      return 0;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->capacity;
}

/** changes capacity of the knapsack constraint
 *
 *  @note This method can only be called during problem creation stage (SCIP_STAGE_PROBLEM)
 */
SCIP_RETCODE SCIPchgCapacityKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Longint          capacity            /**< new capacity of knapsack */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("method can only be called during problem creation stage\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->capacity = capacity;

   return SCIP_OKAY;
}

/** gets the number of items in the knapsack constraint */
int SCIPgetNVarsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets the array of variables in the knapsack constraint; the user must not modify this array! */
SCIP_VAR** SCIPgetVarsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the array of weights in the knapsack constraint; the user must not modify this array! */
SCIP_Longint* SCIPgetWeightsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->weights;
}

/** gets the dual solution of the knapsack constraint in the current LP */
SCIP_Real SCIPgetDualsolKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the knapsack constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack constraint\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given knapsack constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a knapsack\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}
