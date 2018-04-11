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

/**@file   heur_alns.c
 * @brief  Adaptive large neighborhood search heuristic that orchestrates popular LNS heuristics
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif
#include "scip/heur_alns.h"
#include "scipdefplugins.h"

#define HEUR_NAME             "alns"
#define HEUR_DESC             "Large neighborhood search heuristic that orchestrates the popular neighborhoods Local Branching, RINS, RENS, DINS etc."
#define HEUR_DISPCHAR         'L'
#define HEUR_PRIORITY         -1100500
#define HEUR_FREQ             20
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define NNEIGHBORHOODS 8

/*
 * limit parameters for sub-SCIPs
 */
#define DEFAULT_NODESQUOT        0.1
#define DEFAULT_NODESOFFSET      500LL
#define DEFAULT_NSOLSLIM         3
#define DEFAULT_MINNODES         50LL
#define DEFAULT_MAXNODES         5000LL
#define DEFAULT_WAITINGNODES     25LL  /**< number of nodes since last incumbent solution that the heuristic should wait */
#define DEFAULT_TARGETNODEFACTOR 1.5
#define DEFAULT_STALLNODEFACTOR  0.25
#define LPLIMFAC                 4.0

/*
 * parameters for the minimum improvement
 */
#define DEFAULT_MINIMPROVELOW    0.0001
#define DEFAULT_MINIMPROVEHIGH   0.1
#define MINIMPROVEFAC            1.5
#define DEFAULT_STARTMINIMPROVE  0.05
#define DEFAULT_ADJUSTMINIMPROVE TRUE

/*
 * bandit algorithm parameters
 */
#define DEFAULT_BESTSOLWEIGHT  1
#define DEFAULT_BANDITALGO     'u'  /**< the default bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy */
#define DEFAULT_GAMMA          0.2  /**< default weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3 */
#define DEFAULT_BETA           0.0  /**< default reward offset between 0 and 1 at every observation for exp3 */
#define DEFAULT_REWARDCONTROL  0.8  /**< reward control to increase the weight of the simple solution indicator and decrease the weight of the closed gap reward */
#define DEFAULT_SCALEBYEFFORT  TRUE /**< should the reward be scaled by the effort? */
#define DEFAULT_EPS            0.5  /**< increase exploration in epsilon-greedy bandit algorithm */
#define DEFAULT_RESETWEIGHTS   TRUE /**< should the bandit algorithms be reset when a new problem is read? */
#define DEFAULT_SUBSCIPRANDSEEDS FALSE /**< should random seeds of sub-SCIPs be altered to increase diversification? */
#define DEFAULT_ALPHA          0.2  /**< parameter to increase the confidence width in UCB */
#define DEFAULT_REWARDBASELINE 0.5  /**< the reward baseline to separate successful and failed calls */
#define DEFAULT_FIXTOL         0.1  /**< tolerance by which the fixing rate may be missed/exceeded without generic (unfixing) */
#define DEFAULT_USELOCALREDCOST FALSE /**< should local reduced costs be used for generic (un)fixing? */

/*
 * parameters to control variable fixing
 */
#define DEFAULT_USEREDCOST       TRUE  /**< should reduced cost scores be used for variable priorization? */
#define DEFAULT_USEPSCOST        TRUE  /**< should pseudo cost scores be used for variable priorization? */
#define DEFAULT_USEDISTANCES     TRUE  /**< should distances from fixed variables be used for variable priorization */
#define DEFAULT_DOMOREFIXINGS    TRUE  /**< should the ALNS heuristic do more fixings by itself based on variable prioritization
                                         *  until the target fixing rate is reached? */
#define DEFAULT_ADJUSTFIXINGRATE TRUE  /**< should the heuristic adjust the target fixing rate based on the success? */
#define FIXINGRATE_DECAY         0.75  /**< geometric decay for fixing rate adjustments */
#define FIXINGRATE_STARTINC      0.2   /**< initial increment value for fixing rate */
#define DEFAULT_USESUBSCIPHEURS  FALSE /**< should the heuristic activate other sub-SCIP heuristics during its search?  */
#define DEFAULT_COPYCUTS         FALSE /**< should cutting planes be copied to the sub-SCIP? */
#define DEFAULT_REWARDFILENAME   "-"   /**< file name to store all rewards and the selection of the bandit */

/* individual random seeds */
#define DEFAULT_SEED 113
#define MUTATIONSEED 121
#define CROSSOVERSEED 321

/* individual neighborhood parameters */
#define DEFAULT_MINFIXINGRATE_RENS 0.3
#define DEFAULT_MAXFIXINGRATE_RENS 0.7
#define DEFAULT_ACTIVE_RENS TRUE
#define DEFAULT_PRIORITY_RENS 1.0

#define DEFAULT_MINFIXINGRATE_RINS 0.2
#define DEFAULT_MAXFIXINGRATE_RINS 0.6
#define DEFAULT_ACTIVE_RINS TRUE
#define DEFAULT_PRIORITY_RINS 1.0

#define DEFAULT_MINFIXINGRATE_MUTATION 0.4
#define DEFAULT_MAXFIXINGRATE_MUTATION 0.9
#define DEFAULT_ACTIVE_MUTATION TRUE
#define DEFAULT_PRIORITY_MUTATION 1.0

#define DEFAULT_MINFIXINGRATE_LOCALBRANCHING 0.0
#define DEFAULT_MAXFIXINGRATE_LOCALBRANCHING 0.9
#define DEFAULT_ACTIVE_LOCALBRANCHING TRUE
#define DEFAULT_PRIORITY_LOCALBRANCHING 1.0

#define DEFAULT_MINFIXINGRATE_PROXIMITY 0.0
#define DEFAULT_MAXFIXINGRATE_PROXIMITY 0.9
#define DEFAULT_ACTIVE_PROXIMITY TRUE
#define DEFAULT_PRIORITY_PROXIMITY 1.0

#define DEFAULT_MINFIXINGRATE_CROSSOVER 0.4
#define DEFAULT_MAXFIXINGRATE_CROSSOVER 0.9
#define DEFAULT_ACTIVE_CROSSOVER TRUE
#define DEFAULT_PRIORITY_CROSSOVER 1.0

#define DEFAULT_MINFIXINGRATE_ZEROOBJECTIVE 0.0
#define DEFAULT_MAXFIXINGRATE_ZEROOBJECTIVE 0.9
#define DEFAULT_ACTIVE_ZEROOBJECTIVE TRUE
#define DEFAULT_PRIORITY_ZEROOBJECTIVE 1.0

#define DEFAULT_MINFIXINGRATE_DINS 0.1
#define DEFAULT_MAXFIXINGRATE_DINS 0.5
#define DEFAULT_ACTIVE_DINS TRUE
#define DEFAULT_PRIORITY_DINS 1.0

#define DEFAULT_NSOLS_CROSSOVER 2 /**< parameter for the number of solutions that crossover should combine */
#define DEFAULT_NPOOLSOLS_DINS  5 /**< number of pool solutions where binary solution values must agree */

/* event handler properties */
#define EVENTHDLR_NAME         "Alns"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"
#define SCIP_EVENTTYPE_ALNS (SCIP_EVENTTYPE_LPSOLVED | SCIP_EVENTTYPE_SOLFOUND | SCIP_EVENTTYPE_BESTSOLFOUND)

/* properties of the ALNS neighborhood statistics table */
#define TABLE_NAME_NEIGHBORHOOD                  "neighborhood"
#define TABLE_DESC_NEIGHBORHOOD                  "ALNS neighborhood statistics"
#define TABLE_POSITION_NEIGHBORHOOD              12500                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_NEIGHBORHOOD        SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

/*
 * Data structures
 */

/*
 * additional neighborhood data structures
 */


typedef struct data_crossover DATA_CROSSOVER; /**< crossover neighborhood data structure */

typedef struct data_mutation DATA_MUTATION; /**< mutation neighborhood data structure */

typedef struct data_dins DATA_DINS; /**< dins neighborhood data structure */

typedef struct NH_FixingRate NH_FIXINGRATE; /** fixing rate data structure */

typedef struct NH_Stats NH_STATS; /**< neighborhood statistics data structure */

typedef struct Nh NH;             /**< neighborhood data structure */


/*
 * variable priorization data structure for sorting
 */
typedef struct VarPrio VARPRIO;

/** callback to collect variable fixings of neighborhood */
 #define DECL_VARFIXINGS(x) SCIP_RETCODE x ( \
   SCIP*                 scip,               /**< SCIP data structure */                     \
   NH*                   neighborhood,       /**< ALNS neighborhood data structure */         \
   SCIP_VAR**            varbuf,             /**< buffer array to collect variables to fix */\
   SCIP_Real*            valbuf,             /**< buffer array to collect fixing values */   \
   int*                  nfixings,           /**< pointer to store the number of fixings */  \
   SCIP_RESULT*          result              /**< result pointer */                          \
   )

/** callback for subproblem changes other than variable fixings
 *
 *  this callback can be used to further modify the subproblem by changes other than variable fixings.
 *  Typical modifications include restrictions of variable domains, the formulation of additional constraints,
 *  or changed objective coefficients.
 *
 *  The callback should set the \p success pointer to indicate whether it was successful with its modifications or not.
 */
#define DECL_CHANGESUBSCIP(x) SCIP_RETCODE x (  \
   SCIP*                 sourcescip,         /**< source SCIP data structure */\
   SCIP*                 targetscip,         /**< target SCIP data structure */\
   SCIP_VAR**            subvars,            /**< array of targetscip variables in the same order as the source SCIP variables */\
   int*                  ndomchgs,           /**< pointer to store the number of performed domain changes */\
   int*                  nchgobjs,           /**< pointer to store the number of changed objective coefficients */ \
   int*                  naddedconss,        /**< pointer to store the number of additional constraints */\
   SCIP_Bool*            success             /**< pointer to store if the sub-MIP was successfully adjusted */\
   )

/** optional initialization callback for neighborhoods when a new problem is read */
#define DECL_NHINIT(x) SCIP_RETCODE x (                                          \
   SCIP*                 scip,               /**< SCIP data structure */         \
   NH*                   neighborhood        /**< neighborhood data structure */ \
   )

/** deinitialization callback for neighborhoods when exiting a problem */
#define DECL_NHEXIT(x) SCIP_RETCODE x ( \
   SCIP*                 scip,               /**< SCIP data structure */         \
   NH*                   neighborhood        /**< neighborhood data structure */ \
   )

/** deinitialization callback for neighborhoods before SCIP is freed */
#define DECL_NHFREE(x) SCIP_RETCODE x (      \
   SCIP*                 scip,               /**< SCIP data structure */         \
   NH*                   neighborhood        /**< neighborhood data structure */ \
   )

/** callback function to return a feasible reference solution for further fixings
 *
 *  The reference solution should be stored in the \p solptr.
 *  The \p result pointer can be used to indicate either
 *
 *  - SCIP_SUCCESS or
 *  - SCIP_DIDNOTFIND
 */
#define DECL_NHREFSOL(x) SCIP_RETCODE x (                                       \
   SCIP*                 scip,               /**< SCIP data structure */  \
   NH*                   neighborhood,       /**< neighborhood data structure */ \
   SCIP_SOL**            solptr,             /**< pointer to store the reference solution */ \
   SCIP_RESULT*          result              /**< pointer to indicate the callback success whether a reference solution is available */ \
   )

/** callback function to deactivate neighborhoods on problems where they are irrelevant */
#define DECL_NHDEACTIVATE(x) SCIP_RETCODE x (\
   SCIP*                 scip,               /**< SCIP data structure */  \
   SCIP_Bool*            deactivate          /**< pointer to store whether the neighborhood should be deactivated (TRUE) for an instance */ \
   )

/** sub-SCIP status code enumerator */
enum HistIndex
{
   HIDX_OPT              = 0,                /**< sub-SCIP was solved to optimality  */
   HIDX_USR              = 1,                /**< sub-SCIP was user interrupted */
   HIDX_NODELIM          = 2,                /**< sub-SCIP reached the node limit */
   HIDX_STALLNODE        = 3,                /**< sub-SCIP reached the stall node limit */
   HIDX_INFEAS           = 4,                /**< sub-SCIP was infeasible */
   HIDX_SOLLIM           = 5,                /**< sub-SCIP reached the solution limit */
   HIDX_OTHER            = 6                 /**< sub-SCIP reached none of the above codes */
};
typedef enum HistIndex HISTINDEX;
#define NHISTENTRIES 7


/** statistics for a neighborhood */
struct NH_Stats
{
   SCIP_CLOCK*           setupclock;         /**< clock for sub-SCIP setup time */
   SCIP_CLOCK*           submipclock;        /**< clock for the sub-SCIP solve */
   SCIP_Longint          usednodes;          /**< total number of used nodes */
   SCIP_Real             oldupperbound;      /**< upper bound before the sub-SCIP started */
   int                   nruns;              /**< number of runs of a neighborhood */
   int                   nrunsbestsol;       /**< number of runs that produced a new incumbent */
   SCIP_Longint          nsolsfound;         /**< the total number of solutions found */
   SCIP_Longint          nbestsolsfound;     /**< the total number of improving solutions found */
   int                   nfixings;           /**< the number of fixings in one run */
   int                   statushist[NHISTENTRIES]; /**< array to count sub-SCIP statuses */
};


/** fixing rate data structure to control the amount of target fixings of a neighborhood */
struct NH_FixingRate
{
   SCIP_Real             minfixingrate;      /**< the minimum fixing rate */
   SCIP_Real             targetfixingrate;   /**< the current target fixing rate */
   SCIP_Real             increment;          /**< the current increment by which the target fixing rate is in-/decreased */
   SCIP_Real             maxfixingrate;      /**< the maximum fixing rate */
};

/** neighborhood data structure with callbacks, statistics, fixing rate */
struct Nh
{
   char*                 name;               /**< the name of this neighborhood */
   NH_FIXINGRATE         fixingrate;         /**< fixing rate for this neighborhood */
   NH_STATS              stats;              /**< statistics for this neighborhood */
   DECL_VARFIXINGS       ((*varfixings));    /**< variable fixings callback for this neighborhood */
   DECL_CHANGESUBSCIP    ((*changesubscip)); /**< callback for subproblem changes other than variable fixings */
   DECL_NHINIT           ((*nhinit));        /**< initialization callback when a new problem is read */
   DECL_NHEXIT           ((*nhexit));        /**< deinitialization callback when exiting a problem */
   DECL_NHFREE           ((*nhfree));        /**< deinitialization callback before SCIP is freed */
   DECL_NHREFSOL         ((*nhrefsol));      /**< callback function to return a reference solution for further fixings, or NULL */
   DECL_NHDEACTIVATE     ((*nhdeactivate));  /**< callback function to deactivate neighborhoods on problems where they are irrelevant, or NULL if it is always active */
   SCIP_Bool             active;             /**< is this neighborhood active or not? */
   SCIP_Real             priority;           /**< positive call priority to initialize bandit algorithms */
   union
   {
      DATA_MUTATION*     mutation;           /**< mutation data */
      DATA_CROSSOVER*    crossover;          /**< crossover data */
      DATA_DINS*         dins;               /**< dins data */
   }                     data;               /**< data object for neighborhood specific data */
};

/** mutation neighborhood data structure */
struct data_mutation
{
   SCIP_RANDNUMGEN*      rng;                /**< random number generator */
};

/** crossover neighborhood data structure */
struct data_crossover
{
   int                   nsols;              /**< the number of solutions that crossover should combine */
   SCIP_RANDNUMGEN*      rng;                /**< random number generator to draw from the solution pool */
   SCIP_SOL*             selsol;             /**< best selected solution by crossover as reference point */
};

/** dins neighborhood data structure */
struct data_dins
{
   int                   npoolsols;          /**< number of pool solutions where binary solution values must agree */
};

/** primal heuristic data */
struct SCIP_HeurData
{
   NH**                  neighborhoods;      /**< array of neighborhoods */
   SCIP_BANDIT*          bandit;             /**< bandit algorithm */
   char*                 rewardfilename;     /**< file name to store all rewards and the selection of the bandit */
   FILE*                 rewardfile;         /**< reward file pointer, or NULL */
   SCIP_Longint          nodesoffset;        /**< offset added to the nodes budget */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes in a single sub-SCIP */
   SCIP_Longint          targetnodes;        /**< targeted number of nodes to start a sub-SCIP */
   SCIP_Longint          minnodes;           /**< minimum number of nodes required to start a sub-SCIP */
   SCIP_Longint          usednodes;          /**< total number of nodes already spent in sub-SCIPs */
   SCIP_Longint          waitingnodes;       /**< number of nodes since last incumbent solution that the heuristic should wait */
   SCIP_Real             nodesquot;          /**< fraction of nodes compared to the main SCIP for budget computation */
   SCIP_Real             startminimprove;    /**< initial factor by which ALNS should at least improve the incumbent */
   SCIP_Real             minimprovelow;      /**< lower threshold for the minimal improvement over the incumbent */
   SCIP_Real             minimprovehigh;     /**< upper bound for the minimal improvement over the incumbent */
   SCIP_Real             minimprove;         /**< factor by which ALNS should at least improve the incumbent */
   SCIP_Real             lplimfac;           /**< limit fraction of LPs per node to interrupt sub-SCIP */
   SCIP_Real             exp3_gamma;         /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3 */
   SCIP_Real             exp3_beta;          /**< reward offset between 0 and 1 at every observation for exp3 */
   SCIP_Real             epsgreedy_eps;      /**< increase exploration in epsilon-greedy bandit algorithm */
   SCIP_Real             ucb_alpha;          /**< parameter to increase the confidence width in UCB */
   SCIP_Real             rewardcontrol;      /**< reward control to increase the weight of the simple solution indicator
                                               *  and decrease the weight of the closed gap reward */
   SCIP_Real             targetnodefactor;   /**< factor by which target node number is increased/decreased at every adjustment */
   SCIP_Real             stallnodefactor;    /**< stall node limit as a fraction of total node limit */
   SCIP_Real             rewardbaseline;     /**< the reward baseline to separate successful and failed calls */
   SCIP_Real             fixtol;             /**< tolerance by which the fixing rate may be missed/exceeded without generic (unfixing) */
   int                   nneighborhoods;     /**< number of neighborhoods */
   int                   nactiveneighborhoods;/**< number of active neighborhoods */
   int                   ninitneighborhoods; /**< neighborhoods that were used at least one time */
   int                   nsolslim;           /**< limit on the number of improving solutions in a sub-SCIP call */
   int                   seed;               /**< initial random seed for bandit algorithms and random decisions by neighborhoods */
   int                   currneighborhood;   /**< index of currently selected neighborhood */
   int                   ndelayedcalls;      /**< the number of delayed calls */
   char                  banditalgo;         /**< the bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy */
   SCIP_Bool             useredcost;         /**< should reduced cost scores be used for variable prioritization? */
   SCIP_Bool             usedistances;       /**< should distances from fixed variables be used for variable prioritization */
   SCIP_Bool             usepscost;          /**< should pseudo cost scores be used for variable prioritization? */
   SCIP_Bool             domorefixings;      /**< should the ALNS heuristic do more fixings by itself based on variable prioritization
                                               *  until the target fixing rate is reached? */
   SCIP_Bool             adjustfixingrate;   /**< should the heuristic adjust the target fixing rate based on the success? */
   SCIP_Bool             usesubscipheurs;    /**< should the heuristic activate other sub-SCIP heuristics during its search?  */
   SCIP_Bool             adjustminimprove;   /**< should the factor by which the minimum improvement is bound be dynamically updated? */
   SCIP_Bool             resetweights;       /**< should the bandit algorithms be reset when a new problem is read? */
   SCIP_Bool             subsciprandseeds;   /**< should random seeds of sub-SCIPs be altered to increase diversification? */
   SCIP_Bool             scalebyeffort;      /**< should the reward be scaled by the effort? */
   SCIP_Bool             copycuts;           /**< should cutting planes be copied to the sub-SCIP? */
   SCIP_Bool             uselocalredcost;    /**< should local reduced costs be used for generic (un)fixing? */
};

/** event handler data */
struct SCIP_EventData
{
   SCIP_VAR**            subvars;            /**< the variables of the subproblem */
   SCIP*                 sourcescip;         /**< original SCIP data structure */
   SCIP_HEUR*            heur;               /**< alns heuristic structure */
   SCIP_Longint          nodelimit;          /**< node limit of the run */
   SCIP_Real             lplimfac;           /**< limit fraction of LPs per node to interrupt sub-SCIP */
   NH_STATS*             runstats;           /**< run statistics for the current neighborhood */
   SCIP_Bool             allrewardsmode;     /**< true if solutions should only be checked for reward comparisons */
};

/** represents limits for the sub-SCIP solving process */
struct SolveLimits
{
   SCIP_Longint          nodelimit;          /**< maximum number of solving nodes for the sub-SCIP */
   SCIP_Real             memorylimit;        /**< memory limit for the sub-SCIP */
   SCIP_Real             timelimit;          /**< time limit for the sub-SCIP */
   SCIP_Longint          stallnodes;         /**< maximum number of nodes without (primal) stalling */
};

typedef struct SolveLimits SOLVELIMITS;

/** data structure that can be used for variable prioritization for additional fixings */
struct VarPrio
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_Real*            randscores;         /**< random scores for prioritization */
   int*                  distances;          /**< breadth-first distances from already fixed variables */
   SCIP_Real*            redcostscores;      /**< reduced cost scores for fixing a variable to a reference value */
   SCIP_Real*            pscostscores;       /**< pseudocost scores for fixing a variable to a reference value */
   unsigned int          useredcost:1;       /**< should reduced cost scores be used for variable prioritization? */
   unsigned int          usedistances:1;     /**< should distances from fixed variables be used for variable prioritization */
   unsigned int          usepscost:1;        /**< should pseudo cost scores be used for variable prioritization? */
};

/*
 * Local methods
 */

/** Reset target fixing rate */
static
SCIP_RETCODE resetFixingRate(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_FIXINGRATE*        fixingrate          /**< heuristic fixing rate */
   )
{
   assert(scip != NULL);
   assert(fixingrate != NULL);
   fixingrate->increment = FIXINGRATE_STARTINC;

   /* use the middle between the minimum and the maximum fixing rate */
   fixingrate->targetfixingrate = 0.5 * (fixingrate->minfixingrate + fixingrate->maxfixingrate);

   return SCIP_OKAY;
}

/** reset the currently active neighborhood */
static
void resetCurrentNeighborhood(
   SCIP_HEURDATA*        heurdata
   )
{
   assert(heurdata != NULL);
   heurdata->currneighborhood = -1;
   heurdata->ndelayedcalls = 0;
}

/** update increment for fixing rate */
static
void updateFixingRateIncrement(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->increment *= FIXINGRATE_DECAY;
   fx->increment = MAX(fx->increment, 0.01);
}


/** increase fixing rate
 *
 *  decrease also the rate by which the target fixing rate is adjusted
 */
static
void increaseFixingRate(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->targetfixingrate += fx->increment;
   fx->targetfixingrate = MIN(fx->targetfixingrate, fx->maxfixingrate);
   updateFixingRateIncrement(fx);
}

/** decrease fixing rate
 *
 *  decrease also the rate by which the target fixing rate is adjusted
 */
static
void decreaseFixingRate(
   NH_FIXINGRATE*        fx                  /**< fixing rate */
   )
{
   fx->targetfixingrate -= fx->increment;
   fx->targetfixingrate = MAX(fx->targetfixingrate, fx->minfixingrate);
   updateFixingRateIncrement(fx);
}

/** update fixing rate based on the results of the current run */
static
void updateFixingRate(
   SCIP*                 scip,               /**< SCIP data structure */
   NH*                   neighborhood,       /**< neighborhood */
   SCIP_STATUS           subscipstatus,      /**< status of the sub-SCIP run */
   NH_STATS*             runstats            /**< run statistics for this run */
   )
{
   NH_FIXINGRATE* fx;

   fx = &neighborhood->fixingrate;

   switch (subscipstatus) {
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_INFEASIBLE:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_SOLLIMIT:
         /* decrease the fixing rate (make subproblem harder) */
         decreaseFixingRate(fx);
         break;
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_USERINTERRUPT:
      case SCIP_STATUS_NODELIMIT:
         /* increase the fixing rate (make the subproblem easier) only if no solution was found */
         if( runstats->nbestsolsfound <= 0 )
            increaseFixingRate(fx);
         break;
      /* fall through cases to please lint */
      case SCIP_STATUS_UNKNOWN:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_TIMELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_BESTSOLLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_UNBOUNDED:
      default:
         break;
   }
}

/** increase target node limit */
static
void increaseTargetNodeLimit(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   heurdata->targetnodes = (SCIP_Longint)(heurdata->targetnodes * heurdata->targetnodefactor);
   heurdata->targetnodes = MIN(heurdata->targetnodes, heurdata->maxnodes);
}

/** decrease target node limit */
static
void decreaseTargetNodeLimit(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   heurdata->targetnodes = (SCIP_Longint)(heurdata->targetnodes / heurdata->targetnodefactor);
   heurdata->targetnodes = MAX(heurdata->targetnodes, heurdata->minnodes);
}

/** reset target node limit */
static
void resetTargetNodeLimit(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   heurdata->targetnodes = heurdata->minnodes;
}

/** update target node limit based on the current run results */
static
void updateTargetNodeLimit(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   NH_STATS*             runstats,           /**< statistics of the run */
   SCIP_STATUS           subscipstatus       /**< status of the sub-SCIP run */
   )
{
   switch (subscipstatus) {
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_INFEASIBLE:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_BESTSOLLIMIT:
         /* the subproblem was easy enough -> use smaller limit to speed up SCIP */
         decreaseTargetNodeLimit(heurdata);
         break;
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_NODELIMIT:
         /* the subproblem could be explored more */
         if( runstats->nbestsolsfound == 0 )
            increaseTargetNodeLimit(heurdata);
         break;
      case SCIP_STATUS_USERINTERRUPT:
      case SCIP_STATUS_UNKNOWN:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_TIMELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_UNBOUNDED:
         break;
      default:
         break;
   }
}

/** reset the minimum improvement for the sub-SCIPs */
static
void resetMinimumImprovement(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   assert(heurdata != NULL);
   heurdata->minimprove = heurdata->startminimprove;
}

/** increase minimum improvement for the sub-SCIPs */
static
void increaseMinimumImprovement(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   assert(heurdata != NULL);

   heurdata->minimprove *= MINIMPROVEFAC;
   heurdata->minimprove = MIN(heurdata->minimprove, heurdata->minimprovehigh);
}

/** decrease the minimum improvement for the sub-SCIPs */
static
void decreaseMinimumImprovement(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   assert(heurdata != NULL);

   heurdata->minimprove /= MINIMPROVEFAC;
   SCIPdebugMessage("%.4f", heurdata->minimprovelow);
   heurdata->minimprove = MAX(heurdata->minimprove, heurdata->minimprovelow);
}

/** update the minimum improvement based on the status of the sub-SCIP */
static
void updateMinimumImprovement(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_STATUS           subscipstatus,      /**< status of the sub-SCIP run */
   NH_STATS*             runstats            /**< run statistics for this run */
   )
{
   assert(heurdata != NULL);

   /* if the sub-SCIP status was infeasible, we rather want to make the sub-SCIP easier
    * with a smaller minimum improvement.
    *
    * If a solution limit was reached, we may, set it higher.
    */
   switch (subscipstatus) {
      case SCIP_STATUS_INFEASIBLE:
      case SCIP_STATUS_INFORUNBD:
         /* subproblem was infeasible, probably due to the minimum improvement -> decrease minimum improvement */
         decreaseMinimumImprovement(heurdata);

         break;
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_BESTSOLLIMIT:
      case SCIP_STATUS_OPTIMAL:
         /* subproblem could be optimally solved -> try higher minimum improvement */
         increaseMinimumImprovement(heurdata);
         break;
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_USERINTERRUPT:
         /* subproblem was too hard, decrease minimum improvement */
         if( runstats->nbestsolsfound <= 0 )
            decreaseMinimumImprovement(heurdata);
         break;
      case SCIP_STATUS_UNKNOWN:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_TIMELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_UNBOUNDED:
      default:
         break;
   }
}

/** Reset neighborhood statistics */
static
SCIP_RETCODE neighborhoodStatsReset(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             stats               /**< neighborhood statistics */
   )
{
   assert(scip != NULL);
   assert(stats != NULL);

   stats->nbestsolsfound = 0;
   stats->nruns = 0;
   stats->nrunsbestsol = 0;
   stats->nsolsfound = 0;
   stats->usednodes = 0L;
   stats->nfixings = 0L;

   BMSclearMemoryArray(stats->statushist, NHISTENTRIES);

   SCIP_CALL( SCIPresetClock(scip, stats->setupclock) );
   SCIP_CALL( SCIPresetClock(scip, stats->submipclock) );

   return SCIP_OKAY;
}

/** create a neighborhood of the specified name and include it into the ALNS heuristic */
static
SCIP_RETCODE alnsIncludeNeighborhood(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the ALNS heuristic */
   NH**                  neighborhood,       /**< pointer to store the neighborhood */
   const char*           name,               /**< name for this neighborhood */
   SCIP_Real             minfixingrate,      /**< default value for minfixingrate parameter of this neighborhood */
   SCIP_Real             maxfixingrate,      /**< default value for maxfixingrate parameter of this neighborhood */
   SCIP_Bool             active,             /**< default value for active parameter of this neighborhood */
   SCIP_Real             priority,           /**< positive call priority to initialize bandit algorithms */
   DECL_VARFIXINGS       ((*varfixings)),    /**< variable fixing callback for this neighborhood, or NULL */
   DECL_CHANGESUBSCIP    ((*changesubscip)), /**< subscip changes callback for this neighborhood, or NULL */
   DECL_NHINIT           ((*nhinit)),        /**< initialization callback for neighborhood, or NULL */
   DECL_NHEXIT           ((*nhexit)),        /**< deinitialization callback for neighborhood, or NULL */
   DECL_NHFREE           ((*nhfree)),        /**< deinitialization callback before SCIP is freed, or NULL */
   DECL_NHREFSOL         ((*nhrefsol)),      /**< callback function to return a reference solution for further fixings, or NULL */
   DECL_NHDEACTIVATE     ((*nhdeactivate))   /**< callback function to deactivate neighborhoods on problems where they are irrelevant, or NULL if neighborhood is always active */
   )
{
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(neighborhood != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, neighborhood) );
   assert(*neighborhood != NULL);

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*neighborhood)->name, name, strlen(name)+1) );

   SCIP_CALL( SCIPcreateClock(scip, &(*neighborhood)->stats.setupclock) );
   SCIP_CALL( SCIPcreateClock(scip, &(*neighborhood)->stats.submipclock) );

   (*neighborhood)->changesubscip = changesubscip;
   (*neighborhood)->varfixings = varfixings;
   (*neighborhood)->nhinit = nhinit;
   (*neighborhood)->nhexit = nhexit;
   (*neighborhood)->nhfree = nhfree;
   (*neighborhood)->nhrefsol = nhrefsol;
   (*neighborhood)->nhdeactivate = nhdeactivate;

   /* add parameters for this neighborhood */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/alns/%s/minfixingrate", name);
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "minimum fixing rate for this neighborhood",
         &(*neighborhood)->fixingrate.minfixingrate, TRUE, minfixingrate, 0.0, 1.0, NULL, NULL) );
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/alns/%s/maxfixingrate", name);
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "maximum fixing rate for this neighborhood",
         &(*neighborhood)->fixingrate.maxfixingrate, TRUE, maxfixingrate, 0.0, 1.0, NULL, NULL) );
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/alns/%s/active", name);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname, "is this neighborhood active?",
         &(*neighborhood)->active, TRUE, active, NULL, NULL) );
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/alns/%s/priority", name);
   SCIP_CALL( SCIPaddRealParam(scip, paramname, "positive call priority to initialize bandit algorithms",
         &(*neighborhood)->priority, TRUE, priority, 1e-2, 1.0, NULL, NULL) );

   /* add the neighborhood to the ALNS heuristic */
   heurdata->neighborhoods[heurdata->nneighborhoods++] = (*neighborhood);

   return SCIP_OKAY;
}

/** release all data and free neighborhood */
static
SCIP_RETCODE alnsFreeNeighborhood(
   SCIP*                 scip,               /**< SCIP data structure */
   NH**                  neighborhood        /**< pointer to neighborhood that should be freed */
   )
{
   NH* nhptr;
   assert(scip != NULL);
   assert(neighborhood != NULL);

   nhptr = *neighborhood;
   assert(nhptr != NULL);

   BMSfreeMemoryArray(&nhptr->name);

   /* release further, neighborhood specific data structures */
   if( nhptr->nhfree != NULL )
   {
      SCIP_CALL( nhptr->nhfree(scip, nhptr) );
   }

   SCIP_CALL( SCIPfreeClock(scip, &nhptr->stats.setupclock) );
   SCIP_CALL( SCIPfreeClock(scip, &nhptr->stats.submipclock) );

   SCIPfreeBlockMemory(scip, neighborhood);
   *neighborhood = NULL;

   return SCIP_OKAY;
}

/** initialize neighborhood specific data */
static
SCIP_RETCODE neighborhoodInit(
   SCIP*                 scip,               /**< SCIP data structure */
   NH*                   neighborhood        /**< neighborhood to initialize */
   )
{
   assert(scip != NULL);
   assert(neighborhood != NULL);

   /* call the init callback of the neighborhood */
   if( neighborhood->nhinit != NULL )
   {
      SCIP_CALL( neighborhood->nhinit(scip, neighborhood) );
   }

   return SCIP_OKAY;
}

/** deinitialize neighborhood specific data */
static
SCIP_RETCODE neighborhoodExit(
   SCIP*                 scip,               /**< SCIP data structure */
   NH*                   neighborhood        /**< neighborhood to initialize */
   )
{
   assert(scip != NULL);
   assert(neighborhood != NULL);

   if( neighborhood->nhexit != NULL )
   {
      SCIP_CALL( neighborhood->nhexit(scip, neighborhood) );
   }

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE transferSolution(
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_EVENTDATA*       eventdata           /**< event handler data */
   )
{
   SCIP*      sourcescip;         /* original SCIP data structure */
   SCIP_VAR** subvars;            /* the variables of the subproblem */
   SCIP_HEUR* heur;               /* alns heuristic structure */
   SCIP_SOL*  subsol;             /* solution of the subproblem */
   SCIP_VAR** vars;               /* the original problem's variables                */
   int        nvars;
   SCIP_SOL*  newsol;             /* solution to be created for the original problem */
   SCIP_Real* subsolvals;         /* solution values of the subproblem               */
   SCIP_Bool  success;
   NH_STATS*  runstats;
   SCIP_SOL*  oldbestsol;

   assert(subscip != NULL);

   subsol = SCIPgetBestSol(subscip);
   assert(subsol != NULL);

   sourcescip = eventdata->sourcescip;
   subvars = eventdata->subvars;
   heur = eventdata->heur;
   runstats = eventdata->runstats;
   assert(sourcescip != NULL);
   assert(sourcescip != subscip);
   assert(heur != NULL);
   assert(subvars != NULL);
   assert(runstats != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(sourcescip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(sourcescip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(sourcescip, newsol, nvars, vars, subsolvals) );

   oldbestsol = SCIPgetBestSol(sourcescip);

   /* in the special, experimental all rewards mode, the solution is only checked for feasibility
    * but not stored
    */
   if( eventdata->allrewardsmode )
   {
      SCIP_CALL( SCIPcheckSol(sourcescip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );

      if( success )
      {
         runstats->nsolsfound++;
         if( SCIPgetSolTransObj(sourcescip, newsol) < SCIPgetCutoffbound(sourcescip) )
            runstats->nbestsolsfound++;
      }

      SCIP_CALL( SCIPfreeSol(sourcescip, &newsol) );
   }
   else
   {
      /* try to add new solution to scip and free it immediately */
      SCIP_CALL( SCIPtrySolFree(sourcescip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );

      if( success )
      {
         runstats->nsolsfound++;
         if( SCIPgetBestSol(sourcescip) != oldbestsol )
            runstats->nbestsolsfound++;
      }
   }

   SCIPfreeBufferArray(sourcescip, &subsolvals);

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/** execution callback of the event handler
 *
 * transfer new solutions or interrupt the solving process manually
 */
static
SCIP_DECL_EVENTEXEC(eventExecAlns)
{
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_ALNS);
   assert(eventdata != NULL);

   /* treat the different atomic events */
   switch( SCIPeventGetType(event) )
   {
      case SCIP_EVENTTYPE_SOLFOUND:
      case SCIP_EVENTTYPE_BESTSOLFOUND:
         /* try to transfer the solution to the original SCIP */
         SCIP_CALL( transferSolution(scip, eventdata) );
         break;
      case SCIP_EVENTTYPE_LPSOLVED:
         /* interrupt solution process of sub-SCIP */
         if( SCIPgetNLPs(scip) > eventdata->lplimfac * eventdata->nodelimit )
         {
            SCIPdebugMsg(scip, "interrupt after  %" SCIP_LONGINT_FORMAT " LPs\n", SCIPgetNLPs(scip));
            SCIP_CALL( SCIPinterruptSolve(scip) );
         }
         break;
      default:
         break;
   }

   return SCIP_OKAY;
}

/** initialize neighborhood statistics before the next run */
static
void initRunStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             stats               /**< run statistics */
   )
{
   stats->nbestsolsfound = 0;
   stats->nsolsfound = 0;
   stats->usednodes = 0L;
   stats->nfixings = 0;
   stats->oldupperbound = SCIPgetUpperbound(scip);
}

/** update run stats after the sub SCIP was solved */
static
void updateRunStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             stats,              /**< run statistics */
   SCIP*                 subscip             /**< sub-SCIP instance, or NULL */
   )
{
   /* treat an untransformed subscip as if none was created */
   if( subscip != NULL && ! SCIPisTransformed(subscip) )
      subscip = NULL;

   stats->usednodes = subscip != NULL ? SCIPgetNNodes(subscip) : 0L;
}

/** get the histogram index for this status */
static
int getHistIndex(
   SCIP_STATUS           subscipstatus       /**< sub-SCIP status */
   )
{
   switch (subscipstatus)
   {
      case SCIP_STATUS_OPTIMAL:
         return (int)HIDX_OPT;
      case SCIP_STATUS_INFEASIBLE:
         return (int)HIDX_INFEAS;
      case SCIP_STATUS_NODELIMIT:
         return (int)HIDX_NODELIM;
      case SCIP_STATUS_STALLNODELIMIT:
         return (int)HIDX_STALLNODE;
      case SCIP_STATUS_SOLLIMIT:
         return (int)HIDX_SOLLIM;
      case SCIP_STATUS_USERINTERRUPT:
         return (int)HIDX_USR;
      default:
         return (int)HIDX_OTHER;
   } /*lint !e788*/
}

/** print neighborhood statistics */
static
void printNeighborhoodStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   int i;
   int j;
   HISTINDEX statusses[] = {HIDX_OPT, HIDX_INFEAS, HIDX_NODELIM, HIDX_STALLNODE, HIDX_SOLLIM, HIDX_USR, HIDX_OTHER};

   SCIPinfoMessage(scip, file, "Neighborhoods      : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %4s %4s %4s %4s %4s %4s %4s %4s\n",
            "Calls", "SetupTime", "SolveTime", "SolveNodes", "Sols", "Best", "Exp3", "EpsGreedy", "UCB", "TgtFixRate",
            "Opt", "Inf", "Node", "Stal", "Sol", "Usr", "Othr", "Actv");


   /* loop over neighborhoods and fill in statistics */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood;
      SCIP_Real proba;
      SCIP_Real ucb;
      SCIP_Real epsgreedyweight;
      neighborhood = heurdata->neighborhoods[i];
      SCIPinfoMessage(scip, file, "  %-17s:", neighborhood->name);
      SCIPinfoMessage(scip, file, " %10d", neighborhood->stats.nruns);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, neighborhood->stats.setupclock) );
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, neighborhood->stats.submipclock) );
      SCIPinfoMessage(scip, file, " %10" SCIP_LONGINT_FORMAT, neighborhood->stats.usednodes );
      SCIPinfoMessage(scip, file, " %10d", neighborhood->stats.nsolsfound);
      SCIPinfoMessage(scip, file, " %10d", neighborhood->stats.nbestsolsfound);

      proba = 0.0;
      ucb = 1.0;
      epsgreedyweight = -1.0;

      if( heurdata->bandit != NULL && i < heurdata->nactiveneighborhoods )
      {
         switch (heurdata->banditalgo) {
            case 'u':
               ucb = SCIPgetConfidenceBoundUcb(heurdata->bandit, i);
               break;
            case 'g':
               epsgreedyweight = SCIPgetWeightsEpsgreedy(heurdata->bandit)[i];
               break;
            case 'e':
               proba = SCIPgetProbabilityExp3(heurdata->bandit, i);
               break;
            default:
               break;
         }
      }

      SCIPinfoMessage(scip, file, " %10.5f", proba);
      SCIPinfoMessage(scip, file, " %10.5f", epsgreedyweight);
      SCIPinfoMessage(scip, file, " %10.5f", ucb);
      SCIPinfoMessage(scip, file, " %10.3f", neighborhood->fixingrate.targetfixingrate);

      /* loop over status histogram */
      for( j = 0; j < NHISTENTRIES; ++j )
         SCIPinfoMessage(scip, file, " %4d", neighborhood->stats.statushist[statusses[j]]);

      SCIPinfoMessage(scip, file, " %4d", i < heurdata->nactiveneighborhoods ? 1 : 0);
      SCIPinfoMessage(scip, file, "\n");
   }
}

/** update the statistics of the neighborhood based on the sub-SCIP run */
static
void updateNeighborhoodStats(
   SCIP*                 scip,               /**< SCIP data structure */
   NH_STATS*             runstats,           /**< run statistics */
   NH*                   neighborhood,       /**< the selected neighborhood */
   SCIP_STATUS           subscipstatus       /**< status of the sub-SCIP solve */
   )
{  /*lint --e{715}*/
   NH_STATS* stats;
   stats = &neighborhood->stats;

   /* copy run statistics into neighborhood statistics */
   stats->nbestsolsfound += runstats->nbestsolsfound;
   stats->nsolsfound += runstats->nsolsfound;
   stats->usednodes += runstats->usednodes;
   stats->nruns += 1;

   if( runstats->nbestsolsfound > 0 )
      stats->nrunsbestsol += DEFAULT_BESTSOLWEIGHT;
   else if( runstats->nsolsfound > 0 )
      stats->nrunsbestsol++;

   /* update the counter for the subscip status */
   ++stats->statushist[getHistIndex(subscipstatus)];
}

/** sort callback for variable pointers using the ALNS variable prioritization
 *
 *  the variable prioritization works hierarchically as follows. A variable
 *  a has the higher priority over b iff
 *
 *  - variable distances should be used and a has a smaller distance than b
 *  - variable reduced costs should be used and a has a smaller score than b
 *  - variable pseudo costs should be used and a has a smaller score than b
 *  - based on previously assigned random scores
 *
 *  @note: distances are context-based. For fixing more variables,
 *  distances are initialized from the already fixed variables.
 *  For unfixing variables, distances are initialized starting
 *  from the unfixed variables
 */
static
SCIP_DECL_SORTINDCOMP(sortIndCompAlns)
{  /*lint --e{715}*/
   VARPRIO* varprio;
   SCIP* scip;

   varprio = (VARPRIO*)dataptr;
   assert(varprio != NULL);
   assert(varprio->randscores != NULL);

   scip = varprio->scip;
   assert(scip != NULL);

   if( ind1 == ind2 )
      return 0;

   /* priority is on distances, if enabled. The variable which is closer in a breadth-first search sense to
    * the already fixed variables has precedence */
   if( varprio->usedistances )
   {
      int dist1;
      int dist2;

      dist1 = varprio->distances[ind1];
      dist2 = varprio->distances[ind2];

      if( dist1 < 0 )
         dist1 = INT_MAX;

      if( dist2 < 0 )
         dist2 = INT_MAX;

      assert(varprio->distances != NULL);
      if( dist1 < dist2 )
         return -1;
      else if( dist1 > dist2 )
         return 1;
   }

   assert(! varprio->usedistances || varprio->distances[ind1] == varprio->distances[ind2]);

   /* if the indices tie considering reduced costs or distances are disabled -> use reduced cost information instead */
   if( varprio->useredcost )
   {
      assert(varprio->redcostscores != NULL);

      if( SCIPisLT(scip, varprio->redcostscores[ind1], varprio->redcostscores[ind2]) )
         return -1;
      else if( SCIPisGT(scip, varprio->redcostscores[ind1], varprio->redcostscores[ind2]) )
         return 1;
   }

   assert(! varprio->useredcost || SCIPisEQ(scip, varprio->redcostscores[ind1], varprio->redcostscores[ind2]));

   /* use pseudo cost scores if reduced costs are disabled or a tie was found */
   if( varprio->usepscost )
   {
      assert(varprio->pscostscores != NULL);

      /* prefer the variable with smaller pseudocost score */
      if( SCIPisLT(scip, varprio->pscostscores[ind1], varprio->pscostscores[ind2]) )
         return -1;
      else if( SCIPisGT(scip, varprio->pscostscores[ind1], varprio->pscostscores[ind2]) )
         return 1;
   }


   if( varprio->randscores[ind1] < varprio->randscores[ind2] )
      return -1;
   else if( varprio->randscores[ind1] > varprio->randscores[ind2] )
      return 1;

   return ind1 - ind2;
}

/** Compute the reduced cost score for this variable in the reference solution */
static
SCIP_Real getVariableRedcostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable for which the score should be computed */
   SCIP_Real             refsolval,          /**< solution value in reference solution */
   SCIP_Bool             uselocalredcost     /**< should local reduced costs be used for generic (un)fixing? */
   )
{
   SCIP_Real bestbound;
   SCIP_Real redcost;
   SCIP_Real score;
   assert(scip != NULL);
   assert(var != NULL);

   /* prefer column variables */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return SCIPinfinity(scip);

   if( ! uselocalredcost )
   {
      redcost = SCIPvarGetBestRootRedcost(var);

      bestbound = SCIPvarGetBestRootSol(var);

      /* using global reduced costs, the two factors yield a nonnegative score within tolerances */
      assert(SCIPisDualfeasZero(scip, redcost)
         || (SCIPisDualfeasNegative(scip, redcost) && ! SCIPisFeasPositive(scip, refsolval - bestbound))
         || (SCIPisDualfeasPositive(scip, redcost) && ! SCIPisFeasNegative(scip, refsolval - bestbound)));

   }
   else
   {
      /* this can be safely asserted here, since the heuristic would not reach this point, otherwise */
      assert(SCIPhasCurrentNodeLP(scip));
      assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

      redcost = SCIPgetVarRedcost(scip, var);

      bestbound = SCIPvarGetLPSol(var);
   }

   assert(! SCIPisInfinity(scip, REALABS(bestbound)));
   assert(SCIPisDualfeasZero(scip, redcost) || SCIPisFeasIntegral(scip, bestbound));

   score = redcost * (refsolval - bestbound);

   /* max out numerical inaccuracies from global scores */
   if( ! uselocalredcost )
      score = MAX(score, 0.0);

   return score;
}

/** get the pseudo cost score of this variable with respect to the reference solution */
static
SCIP_Real getVariablePscostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable for which the score should be computed */
   SCIP_Real             refsolval           /**< solution value in reference solution */
   )
{
   SCIP_Real rootsolval;
   assert(scip != NULL);
   assert(var != NULL);

   /* variables that aren't LP columns have no pseudocost score */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return 0.0;

   rootsolval = SCIPvarGetRootSol(var);

   /* the score is 0.0 if the values are equal */
   if( SCIPisEQ(scip, rootsolval, refsolval) )
      return 0.0;
   else
      return SCIPgetVarPseudocostVal(scip, var, refsolval - rootsolval);
}

/** add variable and solution value to buffer data structure for variable fixings. The method checks if
 *  the value still lies within the variable bounds. The value stays unfixed otherwise.
 */
static
void tryAdd2variableBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< (source) SCIP variable that should be added to the buffer */
   SCIP_Real             val,                /**< fixing value for this variable */
   SCIP_VAR**            varbuf,             /**< variable buffer to store variables that should be fixed */
   SCIP_Real*            valbuf,             /**< value buffer to store fixing values */
   int*                  nfixings,           /**< pointer to number of fixed buffer variables, will be increased by 1 */
   SCIP_Bool             integer             /**< is this an integer variable? */
   )
{
   assert(SCIPisFeasIntegral(scip, val) || ! SCIPvarIsIntegral(var));
   assert(*nfixings < SCIPgetNVars(scip));

   /* round the value to its nearest integer */
   if( integer )
      val = SCIPfloor(scip, val + 0.5);

   /* only add fixing if it is still valid within the global variable bounds. Invalidity
    * of this solution value may come from a dual reduction that was performed after the solution from which
    * this value originated was found
    */
   if( SCIPvarGetLbGlobal(var) <= val && val <= SCIPvarGetUbGlobal(var) )
   {
      varbuf[*nfixings] = var;
      valbuf[*nfixings] = val;
      ++(*nfixings);
   }
}

/** query neighborhood for a reference solution for further fixings */
static
SCIP_RETCODE neighborhoodGetRefsol(
   SCIP*                 scip,               /**< SCIP data structure */
   NH*                   neighborhood,       /**< ALNS neighborhood data structure */
   SCIP_SOL**            solptr              /**< solution pointer */
   )
{
   assert(solptr != NULL);
   assert(scip != NULL);
   assert(neighborhood != NULL);

   *solptr = NULL;
   if( neighborhood->nhrefsol != NULL )
   {
      SCIP_RESULT result;
      SCIP_CALL( neighborhood->nhrefsol(scip, neighborhood, solptr, &result) );

      if( result == SCIP_DIDNOTFIND )
         *solptr = NULL;
      else
         assert(*solptr != NULL);
   }

   return SCIP_OKAY;
}

/** fix additional variables found in feasible reference solution if the ones that the neighborhood found were not enough
 *
 *  use not always the best solution for the values, but a reference solution provided by the neighborhood itself
 *
 *  @note it may happen that the target fixing rate is not completely reached. This is the case if intermediate,
 *  dual reductions render the solution values of the reference solution infeasible for
 *  the current, global variable bounds.
 */
static
SCIP_RETCODE alnsFixMoreVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the ALNS neighborhood */
   SCIP_SOL*             refsol,             /**< feasible reference solution for more variable fixings */
   SCIP_VAR**            varbuf,             /**< buffer array to store variables to fix */
   SCIP_Real*            valbuf,             /**< buffer array to store fixing values */
   int*                  nfixings,           /**< pointer to store the number of fixings */
   int                   ntargetfixings,     /**< number of required target fixings */
   SCIP_Bool*            success             /**< pointer to store whether the target fixings have been successfully reached */
   )
{
   VARPRIO varprio;
   SCIP_VAR** vars;
   SCIP_Real* redcostscores;
   SCIP_Real* pscostscores;
   SCIP_Real* solvals;
   SCIP_RANDNUMGEN* rng;
   SCIP_VAR** unfixedvars;
   SCIP_Bool* isfixed;
   int* distances;
   int* perm;
   SCIP_Real* randscores;
   int nbinvars;
   int nintvars;
   int nbinintvars;
   int nvars;
   int b;
   int nvarstoadd;
   int nunfixedvars;

   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(success != NULL);
   assert(heurdata != NULL);
   assert(refsol != NULL);

   *success = FALSE;

   /* if the user parameter forbids more fixings, return immediately */
   if( ! heurdata->domorefixings )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nbinintvars = nbinvars + nintvars;

   if( ntargetfixings >= nbinintvars )
      return SCIP_OKAY;

   /* determine the number of required additional fixings */
   nvarstoadd = ntargetfixings - *nfixings;
   if( nvarstoadd == 0 )
      return SCIP_OKAY;

   varprio.usedistances = heurdata->usedistances && (*nfixings >= 1);
   varprio.useredcost = heurdata->useredcost;
   varprio.usepscost = heurdata->usepscost;
   varprio.scip = scip;
   rng = SCIPbanditGetRandnumgen(heurdata->bandit);
   assert(rng != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &randscores, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &distances, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcostscores, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isfixed, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &unfixedvars, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pscostscores, nbinintvars) );

   /* initialize variable graph distances from already fixed variables */
   if( varprio.usedistances )
   {
      SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, NULL, varbuf, *nfixings, distances, INT_MAX, INT_MAX, ntargetfixings) );
   }
   else
   {
      /* initialize all equal distances to make them irrelevant */
      BMSclearMemoryArray(distances, nbinintvars);
   }

   BMSclearMemoryArray(isfixed, nbinintvars);

   /* mark binary and integer variables if they are fixed */
   for( b = 0; b < *nfixings; ++b )
   {
      int probindex;

      assert(varbuf[b] != NULL);
      probindex = SCIPvarGetProbindex(varbuf[b]);
      assert(probindex >= 0);

      if( probindex < nbinintvars )
         isfixed[probindex] = TRUE;
   }


   SCIP_CALL( SCIPgetSolVals(scip, refsol, nbinintvars, vars, solvals) );

   /* assign scores to unfixed every discrete variable of the problem */
   nunfixedvars = 0;
   for( b = 0; b < nbinintvars; ++b )
   {
      SCIP_VAR* var = vars[b];

      /* filter fixed variables */
      if( isfixed[b] )
         continue;

      /* filter variables with a solution value outside its global bounds */
      if( solvals[b] < SCIPvarGetLbGlobal(var) - 0.5 || solvals[b] > SCIPvarGetUbGlobal(var) + 0.5 )
         continue;

      redcostscores[nunfixedvars] = getVariableRedcostScore(scip, var, solvals[b], heurdata->uselocalredcost);
      pscostscores[nunfixedvars] = getVariablePscostScore(scip, var, solvals[b]);

      unfixedvars[nunfixedvars] = var;
      perm[nunfixedvars] = nunfixedvars;
      randscores[nunfixedvars] = SCIPrandomGetReal(rng, 0.0, 1.0);


      /* these assignments are based on the fact that nunfixedvars <= b */
      solvals[nunfixedvars] = solvals[b];
      distances[nunfixedvars] = distances[b];

      SCIPdebugMsg(scip, "Var <%s> scores: dist %3d, red cost %15.9g, pscost %15.9g rand %6.4f\n",
         SCIPvarGetName(var), distances[nunfixedvars], redcostscores[nunfixedvars],
         pscostscores[nunfixedvars], randscores[nunfixedvars]);

      nunfixedvars++;
   }

   /* use selection algorithm (order of the variables does not matter) for quickly completing the fixing */
   varprio.randscores = randscores;
   varprio.distances = distances;
   varprio.redcostscores = redcostscores;
   varprio.pscostscores = pscostscores;

   nvarstoadd = MIN(nunfixedvars, nvarstoadd);

   /* select the first nvarstoadd many variables according to the score */
   if( nvarstoadd < nunfixedvars )
      SCIPselectInd(perm, sortIndCompAlns, &varprio, nvarstoadd, nunfixedvars);

   /* loop over the first elements of the selection defined in permutation. They represent the best variables */
   for( b = 0; b < nvarstoadd; ++b )
   {
      int permindex = perm[b];
      assert(permindex >= 0);
      assert(permindex < nunfixedvars);

      tryAdd2variableBuffer(scip, unfixedvars[permindex], solvals[permindex], varbuf, valbuf, nfixings, TRUE);
   }

   *success = TRUE;

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &pscostscores);
   SCIPfreeBufferArray(scip, &unfixedvars);
   SCIPfreeBufferArray(scip, &isfixed);
   SCIPfreeBufferArray(scip, &solvals);
   SCIPfreeBufferArray(scip, &redcostscores);
   SCIPfreeBufferArray(scip, &distances);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &randscores);

   return SCIP_OKAY;
}

/** create the bandit algorithm for the heuristic depending on the user parameter */
static
SCIP_RETCODE createBandit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_Real*            priorities,         /**< call priorities for active neighborhoods */
   unsigned int          initseed            /**< initial random seed */
   )
{

   switch (heurdata->banditalgo) {
      case 'u':
         SCIP_CALL( SCIPcreateBanditUcb(scip, &heurdata->bandit, priorities,
               heurdata->ucb_alpha, heurdata->nactiveneighborhoods, initseed) );
         break;

      case 'e':
         SCIP_CALL( SCIPcreateBanditExp3(scip, &heurdata->bandit, priorities,
               heurdata->exp3_gamma, heurdata->exp3_beta, heurdata->nactiveneighborhoods, initseed) );
         break;

      case 'g':
         SCIP_CALL( SCIPcreateBanditEpsgreedy(scip, &heurdata->bandit, priorities,
               heurdata->epsgreedy_eps, heurdata->nactiveneighborhoods, initseed) );
         break;

      default:
         SCIPerrorMessage("Unknown bandit parameter %c\n", heurdata->banditalgo);
         return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyAlns)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurAlns(scip) );

   return SCIP_OKAY;
}

/** unfix some of the variables because there are too many fixed
 *
 *  a variable is ideally unfixed if it is close to other unfixed variables
 *  and fixing it has a high reduced cost impact
 */
static
SCIP_RETCODE alnsUnfixVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the ALNS neighborhood */
   SCIP_VAR**            varbuf,             /**< buffer array to store variables to fix */
   SCIP_Real*            valbuf,             /**< buffer array to store fixing values */
   int*                  nfixings,           /**< pointer to store the number of fixings */
   int                   ntargetfixings,     /**< number of required target fixings */
   SCIP_Bool*            success             /**< pointer to store whether the target fixings have been successfully reached */
   )
{
   VARPRIO varprio;
   SCIP_Real* redcostscores;
   SCIP_Real* pscostscores;
   SCIP_Real* randscores;
   SCIP_VAR** unfixedvars;
   SCIP_VAR** varbufcpy;
   SCIP_Real* valbufcpy;
   SCIP_Bool* isfixedvar;
   SCIP_VAR** vars;
   SCIP_RANDNUMGEN* rng;
   int* distances;
   int* fixeddistances;
   int* perm;
   int nvars;
   int i;
   int nbinintvars;
   int nunfixed;

   *success = FALSE;

   nbinintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if( nbinintvars == 0 )
      return SCIP_OKAY;

   assert(*nfixings > 0);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &isfixedvar, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &unfixedvars, nbinintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &distances, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixeddistances, *nfixings) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcostscores, *nfixings) );
   SCIP_CALL( SCIPallocBufferArray(scip, &randscores, *nfixings) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, *nfixings) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pscostscores, *nfixings) );

   SCIP_CALL( SCIPduplicateBufferArray(scip, &varbufcpy, varbuf, *nfixings) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &valbufcpy, valbuf, *nfixings) );

   /*
    * collect the unfixed binary and integer variables
    */
   BMSclearMemoryArray(isfixedvar, nvars);
   /* loop over fixed variables and mark their respective positions as fixed */
   for( i = 0; i < *nfixings; ++i )
   {
      int probindex = SCIPvarGetProbindex(varbuf[i]);

      assert(probindex >= 0);

      isfixedvar[probindex] = TRUE;
   }

   nunfixed = 0;
   vars = SCIPgetVars(scip);
   /* collect unfixed binary and integer variables */
   for( i = 0; i < nbinintvars; ++i )
   {
      if( ! isfixedvar[i] )
         unfixedvars[nunfixed++] = vars[i];
   }

   varprio.usedistances = heurdata->usedistances && nunfixed > 0;

   /* collect distances of all fixed variables from those that are not fixed */
   if( varprio.usedistances )
   {
      SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, NULL, unfixedvars, nunfixed, distances, INT_MAX, INT_MAX, INT_MAX) );

      for( i = 0; i < *nfixings; ++i )
      {
         int probindex = SCIPvarGetProbindex(varbuf[i]);
         if( probindex >= 0 )
            fixeddistances[i] = distances[probindex];
      }
   }
   else
   {
      BMSclearMemoryArray(fixeddistances, *nfixings);
   }

   /* collect reduced cost scores of the fixings and assign random scores */
   rng = SCIPbanditGetRandnumgen(heurdata->bandit);
   for( i = 0; i < *nfixings; ++i )
   {
      SCIP_VAR* fixedvar = varbuf[i];
      SCIP_Real fixval = valbuf[i];

      /* use negative reduced cost and pseudo cost scores to prefer variable fixings with small score */
      redcostscores[i] = - getVariableRedcostScore(scip, fixedvar, fixval, heurdata->uselocalredcost);
      pscostscores[i] = - getVariablePscostScore(scip, fixedvar, fixval);
      randscores[i] = SCIPrandomGetReal(rng, 0.0, 1.0);
      perm[i] = i;

      SCIPdebugMsg(scip, "Var <%s> scores: dist %3d, red cost %15.9g, pscost %15.9g rand %6.4f\n",
            SCIPvarGetName(fixedvar), fixeddistances[i], redcostscores[i], pscostscores[i], randscores[i]);
   }

   varprio.distances = fixeddistances;
   varprio.randscores = randscores;
   varprio.redcostscores = redcostscores;
   varprio.pscostscores = pscostscores;
   varprio.useredcost = heurdata->useredcost;
   varprio.usepscost = heurdata->usepscost;
   varprio.scip = scip;

   /* scores are assigned in such a way that variables with a smaller score should be fixed last */
   SCIPselectDownInd(perm, sortIndCompAlns, &varprio, ntargetfixings, *nfixings);

   /* bring the desired variables to the front of the array */
   for( i = 0; i < ntargetfixings; ++i )
   {
      valbuf[i] = valbufcpy[perm[i]];
      varbuf[i] = varbufcpy[perm[i]];
   }

   *nfixings = ntargetfixings;

   /* free the buffer arrays in reverse order of allocation */
   SCIPfreeBufferArray(scip, &valbufcpy);
   SCIPfreeBufferArray(scip, &varbufcpy);
   SCIPfreeBufferArray(scip, &pscostscores);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &randscores);
   SCIPfreeBufferArray(scip, &redcostscores);
   SCIPfreeBufferArray(scip, &fixeddistances);
   SCIPfreeBufferArray(scip, &distances);
   SCIPfreeBufferArray(scip, &unfixedvars);
   SCIPfreeBufferArray(scip, &isfixedvar);

   *success = TRUE;

   return SCIP_OKAY;
}

/** call variable fixing callback for this neighborhood and orchestrate additional variable fixings, if necessary */
static
SCIP_RETCODE neighborhoodFixVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the ALNS neighborhood */
   NH*                   neighborhood,       /**< neighborhood data structure */
   SCIP_VAR**            varbuf,             /**< buffer array to keep variables that should be fixed */
   SCIP_Real*            valbuf,             /**< buffer array to keep fixing values */
   int*                  nfixings,           /**< pointer to store the number of variable fixings */
   SCIP_RESULT*          result              /**< pointer to store the result of the fixing operation */
   )
{
   int ntargetfixings;

   assert(scip != NULL);
   assert(neighborhood != NULL);
   assert(varbuf != NULL);
   assert(valbuf != NULL);
   assert(nfixings != NULL);
   assert(result != NULL);

   *nfixings = 0;

   *result = SCIP_DIDNOTRUN;

   if( neighborhood->varfixings != NULL )
   {
      SCIP_CALL( neighborhood->varfixings(scip, neighborhood, varbuf, valbuf, nfixings, result) );

      if( *result != SCIP_SUCCESS )
         return SCIP_OKAY;
   }

   assert(neighborhood->varfixings == NULL || *result != SCIP_DIDNOTRUN);

   ntargetfixings = (int)(neighborhood->fixingrate.targetfixingrate * (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip)));
   SCIPdebugMsg(scip, "Neighborhood Fixings/Target: %d / %d\n",*nfixings, ntargetfixings);

   /* if too few fixings, use a strategy to select more variable fixings: randomized, LP graph, ReducedCost based, mix */
   if( (*result == SCIP_SUCCESS || *result == SCIP_DIDNOTRUN) && (*nfixings <= (1.0 - heurdata->fixtol) * ntargetfixings) )
   {
      SCIP_Bool success;
      SCIP_SOL* refsol;

      /* get reference solution from neighborhood */
      SCIP_CALL( neighborhoodGetRefsol(scip, neighborhood, &refsol) );

      /* try to fix more variables based on the reference solution */
      if( refsol != NULL )
      {
         SCIP_CALL( alnsFixMoreVariables(scip, heurdata, refsol, varbuf, valbuf, nfixings, ntargetfixings, &success) );
      }
      else
         success = FALSE;

      if( success )
         *result = SCIP_SUCCESS;
      else if( *result == SCIP_SUCCESS )
         *result = SCIP_DIDNOTFIND;
      else
         *result = SCIP_DIDNOTRUN;

      SCIPdebugMsg(scip, "After additional fixings: %d / %d\n",*nfixings, ntargetfixings);
   }
   else if( (SCIP_Real)(*nfixings) > (1.0 + heurdata->fixtol) * ntargetfixings)
   {
      SCIP_Bool success;

      SCIP_CALL( alnsUnfixVariables(scip, heurdata, varbuf, valbuf, nfixings, ntargetfixings, &success) );

      assert(success);
      *result = SCIP_SUCCESS;
      SCIPdebugMsg(scip, "Unfixed variables, fixed variables remaining: %d\n", ntargetfixings);
   }
   else
   {
      SCIPdebugMsg(scip, "No additional fixings performed\n");
   }

   return SCIP_OKAY;
}

/** change the sub-SCIP by restricting variable domains, changing objective coefficients, or adding constraints */
static
SCIP_RETCODE neighborhoodChangeSubscip(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   NH*                   neighborhood,       /**< neighborhood */
   SCIP_VAR**            targetvars,         /**< array of target SCIP variables aligned with source SCIP variables */
   int*                  ndomchgs,           /**< pointer to store the number of variable domain changes */
   int*                  nchgobjs,           /**< pointer to store the number of changed objective coefficients */
   int*                  naddedconss,        /**< pointer to store the number of added constraints */
   SCIP_Bool*            success             /**< pointer to store whether the sub-SCIP has been successfully modified */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(neighborhood != NULL);
   assert(targetvars != NULL);
   assert(ndomchgs != NULL);
   assert(nchgobjs != NULL);
   assert(naddedconss != NULL);
   assert(success != NULL);

   *success = FALSE;
   *ndomchgs = 0;
   *nchgobjs = 0;
   *naddedconss = 0;

   /* call the change sub-SCIP callback of the neighborhood */
   if( neighborhood->changesubscip != NULL )
   {
      SCIP_CALL( neighborhood->changesubscip(sourcescip, targetscip, targetvars, ndomchgs, nchgobjs, naddedconss, success) );
   }
   else
   {
      *success = TRUE;
   }

   return SCIP_OKAY;
}

/** set sub-SCIP solving limits */
static
SCIP_RETCODE setLimits(
   SCIP*                 subscip,            /**< SCIP data structure */
   SOLVELIMITS*          solvelimits         /**< pointer to solving limits data structure */
   )
{
   assert(subscip != NULL);
   assert(solvelimits != NULL);

   assert(solvelimits->nodelimit >= solvelimits->stallnodes);

   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", solvelimits->nodelimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", solvelimits->stallnodes));
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", solvelimits->timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", solvelimits->memorylimit) );

   return SCIP_OKAY;
}

/** determine limits for a sub-SCIP */
static
SCIP_RETCODE determineLimits(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< this heuristic */
   SOLVELIMITS*          solvelimits,        /**< pointer to solving limits data structure */
   SCIP_Bool*            runagain            /**< can we solve another sub-SCIP with these limits */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real initfactor;
   SCIP_Real nodesquot;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(solvelimits != NULL);
   assert(runagain != NULL);

   heurdata = SCIPheurGetData(heur);

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &solvelimits->timelimit) );
   if( ! SCIPisInfinity(scip, solvelimits->timelimit) )
      solvelimits->timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &solvelimits->memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( ! SCIPisInfinity(scip, solvelimits->memorylimit) )
   {
      solvelimits->memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      solvelimits->memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( solvelimits->timelimit <= 0.0 || solvelimits->memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      *runagain = FALSE;

   nodesquot = heurdata->nodesquot;
   nodesquot *= (SCIPheurGetNBestSolsFound(heur) + 1.0)/(SCIPheurGetNCalls(heur) + 1.0);

   /* calculate the search node limit of the heuristic  */
   solvelimits->nodelimit = (SCIP_Longint)(nodesquot * SCIPgetNNodes(scip));
   solvelimits->nodelimit += heurdata->nodesoffset;
   solvelimits->nodelimit -= heurdata->usednodes;
   solvelimits->nodelimit -= 100 * SCIPheurGetNCalls(heur);
   solvelimits->nodelimit = MIN(heurdata->maxnodes, solvelimits->nodelimit);

   /* use a smaller budget if not all neighborhoods have been initialized yet */
   assert(heurdata->ninitneighborhoods >= 0);
   initfactor = (heurdata->nactiveneighborhoods - heurdata->ninitneighborhoods + 1.0) / (heurdata->nactiveneighborhoods + 1.0);
   solvelimits->nodelimit = (SCIP_Longint)(solvelimits->nodelimit * initfactor);
   solvelimits->stallnodes = (SCIP_Longint)(solvelimits->nodelimit * heurdata->stallnodefactor);

   /* check whether we have enough nodes left to call subproblem solving */
   if( solvelimits->nodelimit < heurdata->targetnodes )
      *runagain = FALSE;

   return SCIP_OKAY;
}

/** return the bandit algorithm that should be used */
static
SCIP_BANDIT* getBandit(
   SCIP_HEURDATA*        heurdata            /**< heuristic data of the ALNS neighborhood */
   )
{
   assert(heurdata != NULL);
   return heurdata->bandit;
}

/** select a neighborhood depending on the selected bandit algorithm */
static
SCIP_RETCODE selectNeighborhood(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the ALNS neighborhood */
   int*                  neighborhoodidx     /**< pointer to store the selected neighborhood index */
   )
{
   SCIP_BANDIT* bandit;
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(neighborhoodidx != NULL);

   *neighborhoodidx = -1;

   bandit = getBandit(heurdata);

   SCIP_CALL( SCIPbanditSelect(bandit, neighborhoodidx) );
   assert(*neighborhoodidx >= 0);

   return SCIP_OKAY;
}

/** Calculate reward based on the selected reward measure */
static
SCIP_RETCODE getReward(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the ALNS neighborhood */
   NH_STATS*             runstats,           /**< run statistics */
   SCIP_Real*            rewardptr           /**< pointer to store the computed reward */
   )
{
   SCIP_Real reward = 0.0;
   SCIP_Real effort;

   assert(runstats->usednodes >= 0);
   assert(runstats->nfixings >= 0);

   /* just add one node to avoid division by zero */
   effort = runstats->usednodes / (SCIP_Real)(heurdata->targetnodes + 1.0);

   /* assume that every fixed variable linearly reduces the subproblem complexity */
   effort = (1.0 - (runstats->nfixings / ((SCIP_Real)SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip)))) * effort;
   assert(rewardptr != NULL);

   /* a positive reward is only assigned if a new incumbent solution was found */
   if( runstats->nbestsolsfound > 0 )
   {
      SCIP_Real bestsolreward;
      SCIP_Real closedgapreward;
      SCIP_Real rewardcontrol = heurdata->rewardcontrol;

      SCIP_Real lb;
      SCIP_Real ub;

      /* the indicator function is simply 1.0 */
      bestsolreward = 1.0;

      ub = SCIPgetUpperbound(scip);
      lb = SCIPgetLowerbound(scip);

      /* compute the closed gap reward */
      if( SCIPisEQ(scip, ub, lb) )
         closedgapreward = 1.0;
      else if( SCIPisInfinity(scip, runstats->oldupperbound) )
         closedgapreward = 1.0;
      else
      {
         closedgapreward = (runstats->oldupperbound - ub) / (runstats->oldupperbound - lb);
      }

      /* the reward is a convex combination of the best solution reward and the reward for the closed gap */
      reward = rewardcontrol * bestsolreward + (1.0 - rewardcontrol) * closedgapreward;

      /* optionally, scale the reward by the involved effort */
      if( heurdata->scalebyeffort )
      {
         /* reward can be larger than 1.0 if a best solution was found within 0 nodes  */
         reward /= (effort + 1.0);
      }

      /* add the baseline and rescale the reward into the interval [baseline, 1.0] */
      reward = heurdata->rewardbaseline + (1.0 - heurdata->rewardbaseline) * reward;
   }
   else
   {
      /* linearly decrease the reward based on the number of nodes spent */
      SCIP_Real maxeffort = heurdata->targetnodes * heurdata->stallnodefactor;
      SCIP_Real usednodes = runstats->usednodes;

      usednodes *= (1.0 - (runstats->nfixings / ((SCIP_Real)SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip))));

      reward = heurdata->rewardbaseline - (usednodes) * heurdata->rewardbaseline / maxeffort;


      reward = MAX(0.0, reward);
   }

   *rewardptr = reward;
   return SCIP_OKAY;
}

/** update internal bandit algorithm statistics for future draws */
static
SCIP_RETCODE updateBanditAlgorithms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data of the ALNS neighborhood */
   SCIP_Real             reward,             /**< measured reward */
   int                   neighborhoodidx     /**< the neighborhood that was chosen */
   )
{
   SCIP_BANDIT* bandit;
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(neighborhoodidx >= 0);
   assert(neighborhoodidx < heurdata->nactiveneighborhoods);

   bandit = getBandit(heurdata);

   SCIPdebugMsg(scip, "Rewarding bandit algorithm action %d with reward %.2f\n", neighborhoodidx, reward);
   SCIP_CALL( SCIPbanditUpdate(bandit, neighborhoodidx, reward) );

   return SCIP_OKAY;
}

/** set up the sub-SCIP parameters, objective cutoff, and solution limits */
static
SCIP_RETCODE setupSubScip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_VAR**            subvars,            /**< array of sub-SCIP variables in the order of the main SCIP */
   SOLVELIMITS*          solvelimits,        /**< pointer to solving limits data structure */
   SCIP_HEUR*            heur,               /**< this heuristic */
   SCIP_Bool             objchgd             /**< did the objective change between the source and the target SCIP? */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real cutoff;
   SCIP_Real upperbound;

   heurdata = SCIPheurGetData(heur);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console unless we are in debug mode */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

#ifdef ALNS_SUBSCIPOUTPUT
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1) );
   /* enable statistic timing inside sub SCIP */
      SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", TRUE) );
#endif

   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->nsolslim) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   if( ! heurdata->usesubscipheurs )
   {
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );
   }

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && ! SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && ! SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* enable conflict analysis and restrict conflict pool */
   if( ! SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( ! SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handlers; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no decutions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && ! SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
   }

   /* add an objective cutoff */
   if( ! SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
   {
      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      if( ! SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
      {
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip)
                            + heurdata->minimprove * SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
         else
            cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
      }
      cutoff = MIN(upperbound, cutoff);

      if( SCIPisObjIntegral(scip) )
         cutoff = SCIPfloor(scip, cutoff);

      SCIPdebugMsg(scip, "Sub-SCIP cutoff: %15.9" SCIP_REAL_FORMAT " (%15.9" SCIP_REAL_FORMAT " in original space)\n",
         cutoff, SCIPretransformObj(scip, cutoff));

      /* if the objective changed between the source and the target SCIP, encode the cutoff as a constraint */
      if( ! objchgd )
      {
         SCIP_CALL(SCIPsetObjlimit(subscip, cutoff));

         SCIPdebugMsg(scip, "Cutoff added as Objective Limit\n");
      }
      else
      {
         SCIP_CONS* objcons;
         int nvars;
         SCIP_VAR** vars;
         int i;

         vars = SCIPgetVars(scip);
         nvars = SCIPgetNVars(scip);

         SCIP_CALL( SCIPcreateConsLinear(subscip, &objcons, "objbound_of_origscip", 0, NULL, NULL, -SCIPinfinity(subscip), cutoff,
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         for( i = 0; i < nvars; ++i)
         {
            if( ! SCIPisFeasZero(subscip, SCIPvarGetObj(vars[i])) )
            {
               SCIP_CALL( SCIPaddCoefLinear(subscip, objcons, subvars[i], SCIPvarGetObj(vars[i])) );
            }
         }
         SCIP_CALL( SCIPaddCons(subscip, objcons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &objcons) );

         SCIPdebugMsg(scip, "Cutoff added as constraint\n");
      }
   }

   /* set solve limits for sub-SCIP */
   SCIP_CALL( setLimits(subscip, solvelimits) );

   /* change random seed of sub-SCIP */
   if( heurdata->subsciprandseeds )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "randomization/randomseedshift", (int)SCIPheurGetNCalls(heur)) );
   }

   SCIPdebugMsg(scip, "Solve Limits: %lld (%lld) nodes (stall nodes), %.1f sec., %d sols\n",
         solvelimits->nodelimit, solvelimits->stallnodes, solvelimits->timelimit, heurdata->nsolslim);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecAlns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** varbuf;
   SCIP_Real* valbuf;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   NH_STATS runstats[NNEIGHBORHOODS];
   SCIP_STATUS subscipstatus[NNEIGHBORHOODS];
   SCIP* subscip = NULL;

   int nfixings;
   int nvars;
   int neighborhoodidx;
   int ntries;
   SCIP_Bool tryagain;
   NH* neighborhood;
   SOLVELIMITS solvelimits;
   SCIP_Bool success;
   SCIP_Bool run;
   SCIP_Bool allrewardsmode;
   SCIP_Real rewards[NNEIGHBORHOODS];
   int banditidx;
   int i;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( heurdata->nactiveneighborhoods == 0 )
      return SCIP_OKAY;

   /* wait for a sufficient number of nodes since last incumbent solution */
   if( SCIPgetDepth(scip) > 0 && SCIPgetBestSol(scip) != NULL
      && (SCIPgetNNodes(scip) - SCIPsolGetNodenum(SCIPgetBestSol(scip))) < heurdata->waitingnodes )
   {
      SCIPdebugMsg(scip, "Waiting nodes not satisfied\n");
      return SCIP_OKAY;
   }

   run = TRUE;
   /* check if budget allows a run of the next selected neighborhood */
   SCIP_CALL( determineLimits(scip, heur, &solvelimits, &run) );
   SCIPdebugMsg(scip, "Budget check: %" SCIP_LONGINT_FORMAT " (%" SCIP_LONGINT_FORMAT ") %s\n", solvelimits.nodelimit, heurdata->targetnodes, run ? "passed" : "must wait");

   if( ! run )
      return SCIP_OKAY;

   /* delay the heuristic if local reduced costs should be used for generic variable unfixing */
   if( heurdata->uselocalredcost && (nodeinfeasible || ! SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL) )
   {
      *result = SCIP_DELAYED;

      return SCIP_OKAY;
   }



   allrewardsmode = heurdata->rewardfile != NULL;

   /* apply some other rules for a fair all rewards mode; in normal execution mode, neighborhoods are iterated through */
   if( allrewardsmode )
   {
      /* most neighborhoods require an incumbent solution */
      if( SCIPgetNSols(scip) < 2 )
      {
         SCIPdebugMsg(scip, "Not enough solutions for all rewards mode\n");
         return SCIP_OKAY;
      }

      /* if the node is infeasible, or has no LP solution, which is required by some neighborhoods
       * if we are not in all rewards mode, the neighborhoods delay themselves individually
       */
      if( nodeinfeasible || ! SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIPdebugMsg(scip, "Delay ALNS heuristic until a feasible node with optimally solved LP relaxation\n");
         *result = SCIP_DELAYED;
         return SCIP_OKAY;
      }
   }

   /* use the neighborhood that requested a delay or select the next neighborhood to run based on the selected bandit algorithm */
   if( heurdata->currneighborhood >= 0 )
   {
      assert(! allrewardsmode);
      banditidx = heurdata->currneighborhood;
      SCIPdebugMsg(scip, "Select delayed neighborhood %d (was delayed %d times)\n", banditidx, heurdata->ndelayedcalls);
   }
   else
   {
      SCIP_CALL( selectNeighborhood(scip, heurdata, &banditidx) );
      SCIPdebugMsg(scip, "Selected neighborhood %d with bandit algorithm\n", banditidx);
   }

   /* in all rewards mode, we simply loop over all heuristics */
   if( ! allrewardsmode )
      neighborhoodidx = banditidx;
   else
      neighborhoodidx = 0;

   assert(neighborhoodidx >= 0);
   assert(heurdata->nactiveneighborhoods > neighborhoodidx);

   /* allocate memory for variable fixings buffer */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuf, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &valbuf, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* initialize neighborhood statistics for a run */
   ntries = 1;
   do
   {
      SCIP_HASHMAP* varmapf;
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_EVENTDATA eventdata;
      int ndomchgs;
      int nchgobjs;
      int naddedconss;
      int v;
      SCIP_RESULT fixresult;
      tryagain = FALSE;
      neighborhood = heurdata->neighborhoods[neighborhoodidx];
      SCIPdebugMsg(scip, "Running '%s' neighborhood %d\n", neighborhood->name, neighborhoodidx);

      initRunStats(scip, &runstats[neighborhoodidx]);
      rewards[neighborhoodidx] = 0.0;

      subscipstatus[neighborhoodidx] = SCIP_STATUS_UNKNOWN;
      SCIP_CALL( SCIPstartClock(scip, neighborhood->stats.setupclock) );

      /* determine variable fixings and objective coefficients of this neighborhood */
      SCIP_CALL( neighborhoodFixVariables(scip, heurdata, neighborhood, varbuf, valbuf, &nfixings, &fixresult) );

      SCIPdebugMsg(scip, "Fix %d/%d variables\n", nfixings, nvars);

      /* Fixing was not successful, either because the fixing rate was not reached (and no additional variable
       * prioritization was used), or the neighborhood requested a delay, e.g., because no LP relaxation solution exists
       * at the current node
       *
       * The ALNS heuristic keeps a delayed neighborhood active and delays itself.
       */
      if( fixresult != SCIP_SUCCESS )
      {
         SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.setupclock) );

         /* to determine all rewards, we cannot delay neighborhoods */
         if( allrewardsmode )
         {
            if( ntries == heurdata->nactiveneighborhoods )
               break;

            neighborhoodidx = (neighborhoodidx + 1) % heurdata->nactiveneighborhoods;
            ntries++;
            tryagain = TRUE;

            continue;
         }


         /* delay the heuristic along with the selected neighborhood
          *
          * if the neighborhood has been delayed for too many consecutive calls, the delay is treated as a failure */
         if( fixresult == SCIP_DELAYED )
         {

            if( heurdata->ndelayedcalls > (SCIPheurGetFreq(heur) / 4 + 1) )
            {
               resetCurrentNeighborhood(heurdata);

               /* use SCIP_DIDNOTFIND to penalize the neighborhood with a bad reward */
               fixresult = SCIP_DIDNOTFIND;
            }
            else if( heurdata->currneighborhood == -1 )
            {
               heurdata->currneighborhood = neighborhoodidx;
               heurdata->ndelayedcalls = 1;
            }
            else
            {
               heurdata->ndelayedcalls++;
            }
         }

         if( fixresult == SCIP_DIDNOTRUN )
         {
            if( ntries < heurdata->nactiveneighborhoods )
            {
               SCIP_CALL( updateBanditAlgorithms(scip, heurdata, 0.0, neighborhoodidx) );
               SCIP_CALL( selectNeighborhood(scip, heurdata, &neighborhoodidx) );
               ntries++;
               tryagain = TRUE;

               SCIPdebugMsg(scip, "Neighborhood cannot run -> try next neighborhood %d\n", neighborhoodidx);
               continue;
            }
            else
               break;
         }

         assert(fixresult == SCIP_DIDNOTFIND || fixresult == SCIP_DELAYED);
         *result = fixresult;
         break;
      }

      *result = SCIP_DIDNOTFIND;

      neighborhood->stats.nfixings += nfixings;
      runstats[neighborhoodidx].nfixings = nfixings;

      SCIP_CALL( SCIPcreate(&subscip) );
      SCIP_CALL( SCIPhashmapCreate(&varmapf, SCIPblkmem(scip), nvars) );

      /* todo later: run global propagation for this set of fixings */
      SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapf, neighborhood->name, varbuf, valbuf, nfixings, FALSE, heurdata->copycuts, &success, NULL) );

      /* store sub-SCIP variables in array for faster access */
      for( v = 0; v < nvars; ++v )
      {
         subvars[v] = (SCIP_VAR*)SCIPhashmapGetImage(varmapf, (void *)vars[v]);
         assert(subvars[v] != NULL);
      }

      SCIPhashmapFree(&varmapf);

      /* let the neighborhood add additional constraints, or restrict domains */
      SCIP_CALL( neighborhoodChangeSubscip(scip, subscip, neighborhood, subvars, &ndomchgs, &nchgobjs, &naddedconss, &success) );

      if( ! success )
      {
         SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.setupclock) );

         if( ! allrewardsmode || ntries == heurdata->nactiveneighborhoods )
            break;

         neighborhoodidx = (neighborhoodidx + 1) % heurdata->nactiveneighborhoods;
         ntries++;
         tryagain = TRUE;

         SCIP_CALL( SCIPfree(&subscip) );

         continue;
      }

      /* set up sub-SCIP parameters */
      SCIP_CALL( setupSubScip(scip, subscip, subvars, &solvelimits, heur, nchgobjs > 0) );

      /* copy the necessary data into the event data to create new solutions */
      eventdata.nodelimit = solvelimits.nodelimit;
      eventdata.lplimfac = heurdata->lplimfac;
      eventdata.heur = heur;
      eventdata.sourcescip = scip;
      eventdata.subvars = subvars;
      eventdata.runstats = &runstats[neighborhoodidx];
      eventdata.allrewardsmode = allrewardsmode;

      /* include an event handler to transfer solutions into the main SCIP */
      SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAlns, NULL) );

      /* transform the problem before catching the events */
      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_ALNS, eventhdlr, &eventdata, NULL) );

      SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.setupclock) );

      /* todo alternatively: set up sub-SCIP and run presolving */
      /* todo was presolving successful enough regarding fixings? otherwise terminate */

      SCIP_CALL( SCIPstartClock(scip, neighborhood->stats.submipclock) );
      /* run sub-SCIP for the given budget, and collect statistics */
      SCIP_CALL_ABORT( SCIPsolve(subscip) );

#ifdef ALNS_SUBSCIPOUTPUT
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );
#endif

      SCIP_CALL( SCIPstopClock(scip, neighborhood->stats.submipclock) );

      /* update statistics based on the sub-SCIP run results */
      updateRunStats(scip, &runstats[neighborhoodidx], subscip);
      subscipstatus[neighborhoodidx] = SCIPgetStatus(subscip);

      SCIP_CALL( getReward(scip, heurdata, &runstats[neighborhoodidx], &rewards[neighborhoodidx]) );

      /* in all rewards mode, continue with the next neighborhood */
      if( allrewardsmode && ntries < heurdata->nactiveneighborhoods )
      {
         neighborhoodidx = (neighborhoodidx + 1) % heurdata->nactiveneighborhoods;
         ntries++;
         tryagain = TRUE;

         SCIP_CALL( SCIPfree(&subscip) );
      }
   }
   while( tryagain && ! SCIPisStopped(scip) );

   if( subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&subscip) );
   }

   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &valbuf);
   SCIPfreeBufferArray(scip, &varbuf);

   /* update bandit index that may have changed unless we are in all rewards mode */
   if( ! allrewardsmode )
      banditidx = neighborhoodidx;

   if( *result != SCIP_DELAYED )
   {
      /* decrease the number of neighborhoods that have not been initialized */
      if( neighborhood->stats.nruns == 0 )
         --heurdata->ninitneighborhoods;

      heurdata->usednodes += runstats[banditidx].usednodes;

      /* determine the success of this neighborhood, and update the target fixing rate for the next time */
      updateNeighborhoodStats(scip, &runstats[banditidx], heurdata->neighborhoods[banditidx], subscipstatus[banditidx]);

      SCIPdebugMsg(scip, "Status of sub-SCIP run: %d\n", subscipstatus[banditidx]);

      /* adjust the fixing rate for this neighborhood
       * make no adjustments in all rewards mode, because this only affects 1 of 8 heuristics
       */
      if( heurdata->adjustfixingrate && ! allrewardsmode )
      {
         SCIPdebugMsg(scip, "Update fixing rate: %.2f\n", heurdata->neighborhoods[banditidx]->fixingrate.targetfixingrate);
         updateFixingRate(scip, heurdata->neighborhoods[banditidx], subscipstatus[banditidx], &runstats[banditidx]);
         SCIPdebugMsg(scip, "New fixing rate: %.2f\n", heurdata->neighborhoods[banditidx]->fixingrate.targetfixingrate);
      }
      /* similarly, update the minimum improvement for the ALNS heuristic
       * make no adjustments in all rewards mode
       */
      if( heurdata->adjustminimprove && ! allrewardsmode )
      {
         SCIPdebugMsg(scip, "Update Minimum Improvement: %.4f\n", heurdata->minimprove);
         updateMinimumImprovement(heurdata, subscipstatus[banditidx], &runstats[banditidx]);
         SCIPdebugMsg(scip, "--> %.4f\n", heurdata->minimprove);
      }

      /* update the target node limit based on the status of the selected algorithm */
      updateTargetNodeLimit(heurdata, &runstats[banditidx], subscipstatus[banditidx]);

      /* update the bandit algorithms by the measured reward */
      SCIP_CALL( updateBanditAlgorithms(scip, heurdata, rewards[banditidx], banditidx) );

      resetCurrentNeighborhood(heurdata);
   }

   /* write single, measured rewards and the bandit index to the reward file */
   if( allrewardsmode )
   {
      for( i = 0; i < heurdata->nactiveneighborhoods; ++i )
      {
         fprintf(heurdata->rewardfile, "%.4f,", rewards[i]);
      }
      fprintf(heurdata->rewardfile, "%d\n", banditidx);
   }

   return SCIP_OKAY;
}

/** callback to collect variable fixings of RENS */
static
DECL_VARFIXINGS(varFixingsRens)
{  /*lint --e{715}*/
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   int i;
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   *result = SCIP_DELAYED;

   if( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get variable information */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* return if no binary or integer variables are present */
   if( nbinvars + nintvars == 0 )
      return SCIP_OKAY;

   /* loop over binary and integer variables; determine those that should be fixed in the sub-SCIP */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var = vars[i];
      SCIP_Real lpsolval = SCIPgetSolVal(scip, NULL, var);
      assert((i < nbinvars && SCIPvarIsBinary(var)) || (i >= nbinvars && SCIPvarIsIntegral(var)));

      /* fix all binary and integer variables with integer LP solution value */
      if( SCIPisFeasIntegral(scip, lpsolval) )
         tryAdd2variableBuffer(scip, var, lpsolval, varbuf, valbuf, nfixings, TRUE);
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** callback for RENS subproblem changes */
static
DECL_CHANGESUBSCIP(changeSubscipRens)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nintvars;
   int nbinvars;
   int i;

   assert(SCIPhasCurrentNodeLP(sourcescip));
   assert(SCIPgetLPSolstat(sourcescip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* get variable information */
   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* restrict bounds of integer variables with fractional solution value */
   for( i = nbinvars; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var = vars[i];
      SCIP_Real lpsolval = SCIPgetSolVal(sourcescip, NULL, var);

      if( ! SCIPisFeasIntegral(sourcescip, lpsolval) )
      {
         SCIP_Real newlb = SCIPfloor(sourcescip, lpsolval);
         SCIP_Real newub = newlb + 1.0;

         /* only count this as a domain change if the new lower and upper bound are a further restriction */
         if( newlb > SCIPvarGetLbGlobal(subvars[i]) + 0.5 || newub < SCIPvarGetUbGlobal(subvars[i]) - 0.5 )
         {
            SCIP_CALL( SCIPchgVarLbGlobal(targetscip, subvars[i], newlb) );
            SCIP_CALL( SCIPchgVarUbGlobal(targetscip, subvars[i], newub) );
            (*ndomchgs)++;
         }
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/** collect fixings by matching solution values in a collection of solutions for all binary and integer variables,
 *  or for a custom set of variables
 */
static
SCIP_RETCODE fixMatchingSolutionValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sols,               /**< array of 2 or more solutions. It is okay for the array to contain one element
                                               *  equal to NULL to represent the current LP solution */
   int                   nsols,              /**< number of solutions in the array */
   SCIP_VAR**            vars,               /**< variable array for which solution values must agree */
   int                   nvars,              /**< number of variables, or -1 for all binary and integer variables */
   SCIP_VAR**            varbuf,             /**< buffer storage for variable fixings */
   SCIP_Real*            valbuf,             /**< buffer storage for fixing values */
   int*                  nfixings            /**< pointer to store the number of fixings */
   )
{
   int v;
   int nbinintvars;
   SCIP_SOL* firstsol;

   assert(scip != NULL);
   assert(sols != NULL);
   assert(nsols >= 2);
   assert(varbuf != NULL);
   assert(valbuf != NULL);
   assert(nfixings != NULL);
   assert(*nfixings == 0);

   if( nvars == -1 || vars == NULL )
   {
      int nbinvars;
      int nintvars;
      SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );
      nbinintvars = nbinvars + nintvars;
      nvars = nbinintvars;
   }
   firstsol = sols[0];
   assert(nvars > 0);

   /* loop over integer and binary variables and check if their solution values match in all solutions */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real solval;
      SCIP_VAR* var;
      int s;

      var = vars[v];
      assert((v < SCIPgetNBinVars(scip) && SCIPvarIsBinary(var)) || (v >= SCIPgetNBinVars(scip) && SCIPvarIsIntegral(var)));
      solval = SCIPgetSolVal(scip, firstsol, var);

      /* determine if solution values match in all given solutions */
      for( s = 1; s < nsols; ++s )
      {
         SCIP_Real solval2 = SCIPgetSolVal(scip, sols[s], var);
         if( ! SCIPisEQ(scip, solval, solval2) )
            break;
      }

      /* if we did not break early, all solutions agree on the solution value of this variable */
      if( s == nsols )
      {
         tryAdd2variableBuffer(scip, var, solval, varbuf, valbuf, nfixings, TRUE);
      }
   }

   return SCIP_OKAY;
}

/** callback to collect variable fixings of RINS */
static
DECL_VARFIXINGS(varFixingsRins)
{
   /*lint --e{715}*/
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   SCIP_SOL* incumbent;
   SCIP_SOL* sols[2];
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   *result = SCIP_DELAYED;

   if( ! SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   incumbent = SCIPgetBestSol(scip);
   if( incumbent == NULL )
      return SCIP_OKAY;

   if( SCIPsolGetOrigin(incumbent) == SCIP_SOLORIGIN_ORIGINAL )
      return SCIP_OKAY;

   /* get variable information */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* return if no binary or integer variables are present */
   if( nbinvars + nintvars == 0 )
      return SCIP_OKAY;

   sols[0] = NULL;
   sols[1] = incumbent;

   SCIP_CALL( fixMatchingSolutionValues(scip, sols, 2, vars, nbinvars + nintvars, varbuf, valbuf, nfixings) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** initialization callback for crossover when a new problem is read */
static
DECL_NHINIT(nhInitCrossover)
{  /*lint --e{715}*/
   DATA_CROSSOVER* data;

   data = neighborhood->data.crossover;
   assert(data != NULL);

   if( data->rng != NULL )
      SCIPfreeRandom(scip, &data->rng);

   data->selsol = NULL;

   SCIP_CALL( SCIPcreateRandom(scip, &data->rng, CROSSOVERSEED + (unsigned int)SCIPgetNVars(scip)) );

   return SCIP_OKAY;
}

/** deinitialization callback for crossover when exiting a problem */
static
DECL_NHEXIT(nhExitCrossover)
{  /*lint --e{715}*/
   DATA_CROSSOVER* data;
   data = neighborhood->data.crossover;

   assert(neighborhood != NULL);
   assert(data->rng != NULL);

   SCIPfreeRandom(scip, &data->rng);

   return SCIP_OKAY;
}

/** deinitialization callback for crossover before SCIP is freed */
static
DECL_NHFREE(nhFreeCrossover)
{  /*lint --e{715}*/
   assert(neighborhood->data.crossover != NULL);
   SCIPfreeBlockMemory(scip, &neighborhood->data.crossover);

   return SCIP_OKAY;
}

/** callback to collect variable fixings of crossover */
static
DECL_VARFIXINGS(varFixingsCrossover)
{  /*lint --e{715}*/
   DATA_CROSSOVER* data;
   SCIP_RANDNUMGEN* rng;
   SCIP_SOL** sols;
   SCIP_SOL** scipsols;
   int nsols;
   int lastdraw;
   assert(scip != NULL);
   assert(varbuf != NULL);
   assert(nfixings != NULL);
   assert(valbuf != NULL);

   data = neighborhood->data.crossover;

   assert(data != NULL);
   nsols = data->nsols;
   data->selsol = NULL;

   *result = SCIP_DIDNOTRUN;

   /* return if the pool has not enough solutions */
   if( nsols > SCIPgetNSols(scip) )
      return SCIP_OKAY;

   /* return if no binary or integer variables are present */
   if( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   rng = data->rng;
   lastdraw = SCIPgetNSols(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &sols, nsols) );
   scipsols = SCIPgetSols(scip);

   /* draw as many solutions from the pool as required by crossover, biased towards
    * better solutions; therefore, the sorting of the solutions by objective is implicitly used
    */
   while( nsols > 0 )
   {
      /* no need for randomization anymore, exactly nsols many solutions remain for the selection */
      if( lastdraw == nsols )
      {
         int s;

         /* fill the remaining slots 0,...,nsols - 1 by the solutions at the same places */
         for( s = 0; s < nsols; ++s )
            sols[s] = scipsols[s];

         nsols = 0;
      }
      else
      {
         int nextdraw;

         assert(nsols < lastdraw);

         /* draw from the lastdraw - nsols many solutions nsols - 1, ... lastdraw - 1 such that nsols many solution */
         nextdraw = SCIPrandomGetInt(rng, nsols - 1, lastdraw - 1);
         assert(nextdraw >= 0);

         sols[nsols - 1] = scipsols[nextdraw];
         nsols--;
         lastdraw = nextdraw;
      }
   }

   SCIP_CALL( fixMatchingSolutionValues(scip, sols, data->nsols, NULL, -1, varbuf, valbuf, nfixings) );

   /* store best selected solution as reference solution */
   data->selsol = sols[0];
   assert(data->selsol != NULL);

   *result = SCIP_SUCCESS;

   SCIPfreeBufferArray(scip, &sols);

   return SCIP_OKAY;
}

/** callback for crossover reference solution */
static
DECL_NHREFSOL(nhRefsolCrossover)
{ /*lint --e{715}*/
   DATA_CROSSOVER* data;

   data = neighborhood->data.crossover;

   if( data->selsol != NULL )
   {
      *solptr = data->selsol;
      *result = SCIP_SUCCESS;
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}

/** initialization callback for mutation when a new problem is read */
static
DECL_NHINIT(nhInitMutation)
{  /*lint --e{715}*/
   DATA_MUTATION* data;
   assert(scip != NULL);
   assert(neighborhood != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &neighborhood->data.mutation) );

   data = neighborhood->data.mutation;
   assert(data != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &data->rng, MUTATIONSEED + (unsigned int)SCIPgetNVars(scip)) );

   return SCIP_OKAY;
}

/** deinitialization callback for mutation when exiting a problem */
static
DECL_NHEXIT(nhExitMutation)
{  /*lint --e{715}*/
   DATA_MUTATION* data;
   assert(scip != NULL);
   assert(neighborhood != NULL);
   data = neighborhood->data.mutation;
   assert(data != NULL);

   SCIPfreeRandom(scip, &data->rng);

   SCIPfreeBlockMemory(scip, &neighborhood->data.mutation);

   return SCIP_OKAY;
}

/** callback to collect variable fixings of mutation */
static
DECL_VARFIXINGS(varFixingsMutation)
{  /*lint --e{715}*/
   SCIP_RANDNUMGEN* rng;

   SCIP_VAR** vars;
   SCIP_VAR** varscpy;
   int i;
   int nvars;
   int nbinvars;
   int nintvars;
   int nbinintvars;
   int ntargetfixings;
   SCIP_SOL* incumbentsol;
   SCIP_Real targetfixingrate;

   assert(scip != NULL);
   assert(neighborhood != NULL);
   assert(neighborhood->data.mutation != NULL);
   assert(neighborhood->data.mutation->rng != NULL);
   rng = neighborhood->data.mutation->rng;

   *result = SCIP_DIDNOTRUN;

   /* get the problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nbinintvars = nbinvars + nintvars;
   if( nbinintvars == 0 )
      return SCIP_OKAY;

   incumbentsol = SCIPgetBestSol(scip);
   if( incumbentsol == NULL )
      return SCIP_OKAY;

   targetfixingrate = neighborhood->fixingrate.targetfixingrate;
   ntargetfixings = (int)(targetfixingrate * nbinintvars) + 1;

   /* don't continue if number of discrete variables is too small to reach target fixing rate */
   if( nbinintvars <= ntargetfixings )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* copy variables into a buffer array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &varscpy, vars, nbinintvars) );

   /* partially perturb the array until the number of target fixings is reached */
   for( i = 0; *nfixings < ntargetfixings && i < nbinintvars; ++i )
   {
      int randint = SCIPrandomGetInt(rng, i, nbinintvars - 1);
      assert(randint < nbinintvars);

      if( randint > i )
      {
         SCIPswapPointers((void**)&varscpy[i], (void**)&varscpy[randint]);
      }
      /* copy the selected variables and their solution values into the buffer */
      tryAdd2variableBuffer(scip, varscpy[i], SCIPgetSolVal(scip, incumbentsol, varscpy[i]), varbuf, valbuf, nfixings, TRUE);
   }

   assert(i == nbinintvars || *nfixings == ntargetfixings);

   /* Not reaching the number of target fixings means that there is a significant fraction (at least 1 - targetfixingrate)
    * of variables for which the incumbent solution value does not lie within the global bounds anymore. This is a nonsuccess
    * for the neighborhood (additional fixings are not possible), which is okay because the incumbent solution is
    * significantly outdated
    */
   if( *nfixings == ntargetfixings )
      *result = SCIP_SUCCESS;

   /* free the buffer array */
   SCIPfreeBufferArray(scip, &varscpy);

   return SCIP_OKAY;
}

/** add local branching constraint */
static
SCIP_RETCODE addLocalBranchingConstraint(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_VAR**            subvars,            /**< array of sub SCIP variables in same order as source SCIP variables */
   int                   distance,           /**< right hand side of the local branching constraint */
   SCIP_Bool*            success,            /**< pointer to store of a local branching constraint has been successfully added */
   int*                  naddedconss         /**< pointer to increase the number of added constraints */
   )
{
   int nbinvars;
   int i;
   SCIP_SOL* referencesol;
   SCIP_CONS* localbranchcons;
   SCIP_VAR** vars;
   SCIP_Real* consvals;
   SCIP_Real rhs;

   assert(sourcescip != NULL);
   assert(*success == FALSE);

   nbinvars = SCIPgetNBinVars(sourcescip);
   vars = SCIPgetVars(sourcescip);


   if( nbinvars <= 3 )
      return SCIP_OKAY;

   referencesol = SCIPgetBestSol(sourcescip);
   if( referencesol == NULL )
      return SCIP_OKAY;

   rhs = (SCIP_Real)distance;
   rhs = MAX(rhs, 2.0);


   SCIP_CALL( SCIPallocBufferArray(sourcescip, &consvals, nbinvars) );

   /* loop over binary variables and fill the local branching constraint */
   for( i = 0; i < nbinvars; ++i )
   {
      if( SCIPisEQ(sourcescip, SCIPgetSolVal(sourcescip, referencesol, vars[i]), 0.0) )
         consvals[i] = 1.0;
      else
      {
         consvals[i] = -1.0;
         rhs -= 1.0;
      }
   }

   /* create the local branching constraint in the target scip */
   SCIP_CALL( SCIPcreateConsBasicLinear(targetscip, &localbranchcons, "localbranch", nbinvars, subvars, consvals, -SCIPinfinity(sourcescip), rhs) );
   SCIP_CALL( SCIPaddCons(targetscip, localbranchcons) );
   SCIP_CALL( SCIPreleaseCons(targetscip, &localbranchcons) );

   *naddedconss = 1;
   *success = TRUE;

   SCIPfreeBufferArray(sourcescip, &consvals);

   return SCIP_OKAY;
}

/** callback for local branching subproblem changes */
static
DECL_CHANGESUBSCIP(changeSubscipLocalbranching)
{  /*lint --e{715}*/

   SCIP_CALL( addLocalBranchingConstraint(sourcescip, targetscip, subvars, (int)(0.2 * SCIPgetNBinVars(sourcescip)), success, naddedconss) );

   return SCIP_OKAY;
}

/** callback for proximity subproblem changes */
static
DECL_CHANGESUBSCIP(changeSubscipProximity)
{  /*lint --e{715}*/
   SCIP_SOL* referencesol;
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int nvars;
   int i;

   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   if( nbinvars == 0 )
      return SCIP_OKAY;

   referencesol = SCIPgetBestSol(sourcescip);
   if( referencesol == NULL )
      return SCIP_OKAY;

   /* loop over binary variables, set objective coefficients based on reference solution in a local branching fashion */
   for( i = 0; i < nbinvars; ++i )
   {
      SCIP_Real newobj;
      if( SCIPgetSolVal(sourcescip, referencesol, vars[i]) < 0.5 )
         newobj = -1.0;
      else
         newobj = 1.0;
      SCIP_CALL( SCIPchgVarObj(targetscip, subvars[i], newobj) );
   }

   /* loop over the remaining variables and change their objective coefficients to 0 */
   for( ; i < nvars; ++i )
   {
      SCIP_CALL( SCIPchgVarObj(targetscip, subvars[i], 0.0) );
   }

   *nchgobjs = nvars;
   *success = TRUE;

   return SCIP_OKAY;
}

/** callback for zeroobjective subproblem changes */
static
DECL_CHANGESUBSCIP(changeSubscipZeroobjective)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(*success == FALSE);

   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* do not run if no objective variables are present */
   if( SCIPgetNObjVars(sourcescip) == 0 )
      return SCIP_OKAY;

   /* loop over the variables and change their objective coefficients to 0 */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPchgVarObj(targetscip, subvars[i], 0.0) );
   }

   *nchgobjs = nvars;
   *success = TRUE;

   return SCIP_OKAY;
}

/** compute tightened bounds for integer variables depending on how much the LP and the incumbent solution values differ */
static
void computeIntegerVariableBoundsDins(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP_VAR*             var,                /**< the variable for which bounds should be computed */
   SCIP_Real*            lbptr,              /**< pointer to store the lower bound in the DINS sub-SCIP */
   SCIP_Real*            ubptr               /**< pointer to store the upper bound in the DINS sub-SCIP */
   )
{
   SCIP_Real mipsol;
   SCIP_Real lpsol;

   SCIP_Real lbglobal;
   SCIP_Real ubglobal;
   SCIP_SOL* bestsol;

   /* get the bounds for each variable */
   lbglobal = SCIPvarGetLbGlobal(var);
   ubglobal = SCIPvarGetUbGlobal(var);

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);
   /* get the current LP solution for each variable */
   lpsol = SCIPvarGetLPSol(var);

   /* get the current MIP solution for each variable */
   bestsol = SCIPgetBestSol(scip);
   mipsol = SCIPgetSolVal(scip, bestsol, var);


   /* if the solution values differ by 0.5 or more, the variable is rebounded, otherwise it is just copied */
   if( REALABS(lpsol - mipsol) >= 0.5 )
   {
      SCIP_Real range;

      *lbptr = lbglobal;
      *ubptr = ubglobal;

      /* create an equally sized range around lpsol for general integers: bounds are lpsol +- (mipsol-lpsol) */
      range = 2 * lpsol - mipsol;

      if( mipsol >= lpsol )
      {
         range = SCIPfeasCeil(scip, range);
         *lbptr = MAX(*lbptr, range);

         /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
         if( SCIPisFeasEQ(scip, mipsol, *lbptr) )
            *ubptr = *lbptr;
         else
            *ubptr = mipsol;
      }
      else
      {
         range = SCIPfeasFloor(scip, range);
         *ubptr = MIN(*ubptr, range);

         /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
         if( SCIPisFeasEQ(scip, mipsol, *ubptr) )
            *lbptr = *ubptr;
         else
            *lbptr = mipsol;
      }

      /* the global domain of variables might have been reduced since incumbent was found: adjust lb and ub accordingly */
      *lbptr = MAX(*lbptr, lbglobal);
      *ubptr = MIN(*ubptr, ubglobal);
   }
   else
   {
      /* the global domain of variables might have been reduced since incumbent was found: adjust it accordingly */
      *lbptr = MAX(mipsol, lbglobal);
      *ubptr = MIN(mipsol, ubglobal);
   }
}

/** callback to collect variable fixings of DINS */
static
DECL_VARFIXINGS(varFixingsDins)
{
   DATA_DINS* data;
   SCIP_SOL* rootlpsol;
   SCIP_SOL** sols;
   int nsols;
   int nmipsols;
   int nbinvars;
   int nintvars;
   SCIP_VAR** vars;
   int v;

   data = neighborhood->data.dins;
   assert(data != NULL);
   nmipsols = SCIPgetNSols(scip);
   nmipsols = MIN(nmipsols, data->npoolsols);

   *result = SCIP_DELAYED;

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   if( nmipsols == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   if( nbinvars + nintvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcreateSol(scip, &rootlpsol, NULL) );

   /* save root solution LP values in solution */
   for( v = 0; v < nbinvars + nintvars; ++v )
   {
      SCIP_CALL( SCIPsetSolVal(scip, rootlpsol, vars[v], SCIPvarGetRootSol(vars[v])) );
   }

   /* add the node and the root LP solution */
   nsols = nmipsols + 2;

   SCIP_CALL( SCIPallocBufferArray(scip, &sols, nsols) );
   sols[0] = NULL; /* node LP solution */
   sols[1] = rootlpsol;

   /* copy the remaining MIP solutions after the LP solutions */
   BMScopyMemoryArray(&sols[2], SCIPgetSols(scip), nmipsols); /*lint !e866*/

   /* 1. Binary variables are fixed if their values agree in all the solutions */
   if( nbinvars > 0 )
   {
      SCIP_CALL( fixMatchingSolutionValues(scip, sols, nsols, vars, nbinvars, varbuf, valbuf, nfixings) );
   }

   /* 2. Integer variables are fixed if they have a very low distance between the incumbent and the root LP solution */
   for( v = nbinvars; v < nintvars; ++v )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      computeIntegerVariableBoundsDins(scip, vars[v], &lb, &ub);

      if( ub - lb < 0.5 )
      {
         assert(SCIPisFeasIntegral(scip, lb));
         tryAdd2variableBuffer(scip, vars[v], lb, varbuf, valbuf, nfixings, TRUE);
      }
   }

   *result = SCIP_SUCCESS;

   SCIPfreeBufferArray(scip, &sols);

   SCIP_CALL( SCIPfreeSol(scip, &rootlpsol) );

   return SCIP_OKAY;
}

/** callback for DINS subproblem changes */
static
DECL_CHANGESUBSCIP(changeSubscipDins)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nintvars;
   int nbinvars;
   int v;

   SCIP_CALL( SCIPgetVarsData(sourcescip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* 1. loop over integer variables and tighten the bounds */
   for( v = nbinvars; v < nintvars; ++v )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      computeIntegerVariableBoundsDins(sourcescip, vars[v], &lb, &ub);

      SCIP_CALL( SCIPchgVarLbGlobal(targetscip, subvars[v], lb) );
      SCIP_CALL( SCIPchgVarUbGlobal(targetscip, subvars[v], ub) );
      ++(*ndomchgs);
   }

   /* 2. add local branching constraint for binary variables */
   SCIP_CALL( addLocalBranchingConstraint(sourcescip, targetscip, subvars, (int)(0.1 * SCIPgetNBinVars(sourcescip)), success, naddedconss) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** deinitialization callback for DINS before SCIP is freed */
static
DECL_NHFREE(nhFreeDins)
{
   assert(neighborhood->data.dins != NULL);

   SCIPfreeBlockMemory(scip, &neighborhood->data.dins);

   return SCIP_OKAY;
}

/** callback that returns the incumbent solution as a reference point */
static
DECL_NHREFSOL(nhRefsolIncumbent)
{  /*lint --e{715}*/
   assert(scip != NULL);

   if( SCIPgetBestSol(scip) != NULL )
   {
      *result = SCIP_SUCCESS;
      *solptr = SCIPgetBestSol(scip);
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}


/** callback function that deactivates a neighborhood on problems with no discrete variables */
static
DECL_NHDEACTIVATE(nhDeactivateDiscreteVars)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(deactivate != NULL);

   /* deactivate if no discrete variables are present */
   *deactivate = (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0);

   return SCIP_OKAY;
}

/** callback function that deactivates a neighborhood on problems with no binary variables */
static
DECL_NHDEACTIVATE(nhDeactivateBinVars)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(deactivate != NULL);

   /* deactivate if no discrete variables are present */
   *deactivate = (SCIPgetNBinVars(scip) == 0);

   return SCIP_OKAY;
}

/** callback function that deactivates a neighborhood on problems with no objective variables */
static
DECL_NHDEACTIVATE(nhDeactivateObjVars)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(deactivate != NULL);

   /* deactivate if no discrete variables are present */
   *deactivate = (SCIPgetNObjVars(scip) == 0);

   return SCIP_OKAY;
}


/** include all neighborhoods */
static
SCIP_RETCODE includeNeighborhoods(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data of the ALNS heuristic */
   )
{
   NH* rens;
   NH* rins;
   NH* mutation;
   NH* localbranching;
   NH* crossover;
   NH* proximity;
   NH* zeroobjective;
   NH* dins;

   heurdata->nneighborhoods = 0;

   /* include the RENS neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &rens, "rens",
         DEFAULT_MINFIXINGRATE_RENS, DEFAULT_MAXFIXINGRATE_RENS, DEFAULT_ACTIVE_RENS, DEFAULT_PRIORITY_RENS,
         varFixingsRens, changeSubscipRens, NULL, NULL, NULL, NULL, nhDeactivateDiscreteVars) );

   /* include the RINS neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &rins, "rins",
         DEFAULT_MINFIXINGRATE_RINS, DEFAULT_MAXFIXINGRATE_RINS, DEFAULT_ACTIVE_RINS, DEFAULT_PRIORITY_RINS,
         varFixingsRins, NULL, NULL, NULL, NULL, nhRefsolIncumbent, nhDeactivateDiscreteVars) );

   /* include the mutation neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &mutation, "mutation",
         DEFAULT_MINFIXINGRATE_MUTATION, DEFAULT_MAXFIXINGRATE_MUTATION, DEFAULT_ACTIVE_MUTATION, DEFAULT_PRIORITY_MUTATION,
         varFixingsMutation, NULL, nhInitMutation, nhExitMutation, NULL, nhRefsolIncumbent, nhDeactivateDiscreteVars) );

   /* include the local branching neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &localbranching, "localbranching",
         DEFAULT_MINFIXINGRATE_LOCALBRANCHING, DEFAULT_MAXFIXINGRATE_LOCALBRANCHING, DEFAULT_ACTIVE_LOCALBRANCHING, DEFAULT_PRIORITY_LOCALBRANCHING,
         NULL, changeSubscipLocalbranching, NULL, NULL, NULL, nhRefsolIncumbent, nhDeactivateBinVars) );

   /* include the crossover neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &crossover, "crossover",
         DEFAULT_MINFIXINGRATE_CROSSOVER, DEFAULT_MAXFIXINGRATE_CROSSOVER, DEFAULT_ACTIVE_CROSSOVER, DEFAULT_PRIORITY_CROSSOVER,
         varFixingsCrossover, NULL,
         nhInitCrossover, nhExitCrossover, nhFreeCrossover, nhRefsolCrossover, nhDeactivateDiscreteVars) );

   /* allocate data for crossover to include the parameter */
   SCIP_CALL( SCIPallocBlockMemory(scip, &crossover->data.crossover) );
   crossover->data.crossover->rng = NULL;

   /* add crossover neighborhood parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/alns/crossover/nsols", "the number of solutions that crossover should combine",
         &crossover->data.crossover->nsols, TRUE, DEFAULT_NSOLS_CROSSOVER, 2, 10, NULL, NULL) );

   /* include the Proximity neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &proximity, "proximity",
         DEFAULT_MINFIXINGRATE_PROXIMITY, DEFAULT_MAXFIXINGRATE_PROXIMITY, DEFAULT_ACTIVE_PROXIMITY, DEFAULT_PRIORITY_PROXIMITY,
         NULL, changeSubscipProximity, NULL, NULL, NULL, nhRefsolIncumbent, nhDeactivateBinVars) );

   /* include the Zeroobjective neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &zeroobjective, "zeroobjective",
         DEFAULT_MINFIXINGRATE_ZEROOBJECTIVE, DEFAULT_MAXFIXINGRATE_ZEROOBJECTIVE, DEFAULT_ACTIVE_ZEROOBJECTIVE, DEFAULT_PRIORITY_ZEROOBJECTIVE,
         NULL, changeSubscipZeroobjective, NULL, NULL, NULL, nhRefsolIncumbent, nhDeactivateObjVars) );

   /* include the DINS neighborhood */
   SCIP_CALL( alnsIncludeNeighborhood(scip, heurdata, &dins, "dins",
         DEFAULT_MINFIXINGRATE_DINS, DEFAULT_MAXFIXINGRATE_DINS, DEFAULT_ACTIVE_DINS, DEFAULT_PRIORITY_DINS,
         varFixingsDins, changeSubscipDins, NULL, NULL, nhFreeDins, nhRefsolIncumbent, nhDeactivateBinVars) );

   /* allocate data for DINS to include the parameter */
   SCIP_CALL( SCIPallocBlockMemory(scip, &dins->data.dins) );

   /* add DINS neighborhood parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/alns/dins/npoolsols",
         "number of pool solutions where binary solution values must agree",
         &dins->data.dins->npoolsols, TRUE, DEFAULT_NPOOLSOLS_DINS, 1, 100, NULL, NULL) );

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitAlns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* reactivate all neighborhoods if a new problem is read in */
   heurdata->nactiveneighborhoods = heurdata->nneighborhoods;

   /* todo initialize neighborhoods for new problem */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood = heurdata->neighborhoods[i];

      SCIP_CALL( neighborhoodInit(scip, neighborhood) );

      SCIP_CALL( resetFixingRate(scip, &neighborhood->fixingrate) );

      SCIP_CALL( neighborhoodStatsReset(scip, &neighborhood->stats) );
   }

   /* open reward file for reading */
   if( strncasecmp(heurdata->rewardfilename, DEFAULT_REWARDFILENAME, strlen(DEFAULT_REWARDFILENAME)) != 0 )
   {
      heurdata->rewardfile = fopen(heurdata->rewardfilename, "w");

      if( heurdata->rewardfile == NULL )
      {
         SCIPerrorMessage("Error: Could not open reward file <%s>\n", heurdata->rewardfilename);
         return SCIP_FILECREATEERROR;
      }

      SCIPdebugMsg(scip, "Writing reward information to <%s>\n", heurdata->rewardfilename);
   }
   else
      heurdata->rewardfile = NULL;

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolAlns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int i;
   SCIP_Real* priorities;
   unsigned int initseed;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->nactiveneighborhoods = heurdata->nneighborhoods;

   SCIP_CALL( SCIPallocBufferArray(scip, &priorities, heurdata->nactiveneighborhoods) );

   /* init neighborhoods for new problem by resetting their statistics and fixing rate */
   for( i = heurdata->nneighborhoods - 1; i >= 0; --i )
   {
      NH* neighborhood = heurdata->neighborhoods[i];
      SCIP_Bool deactivate;

      SCIP_CALL( neighborhood->nhdeactivate(scip, &deactivate) );

      /* disable inactive neighborhoods */
      if( deactivate || ! neighborhood->active )
      {
         if( heurdata->nactiveneighborhoods - 1 > i )
         {
            assert(heurdata->neighborhoods[heurdata->nactiveneighborhoods - 1]->active);
            SCIPswapPointers((void **)&heurdata->neighborhoods[i], (void **)&heurdata->neighborhoods[heurdata->nactiveneighborhoods - 1]);
         }
         heurdata->nactiveneighborhoods--;
      }
   }

   /* collect neighborhood priorities */
   for( i = 0; i < heurdata->nactiveneighborhoods; ++i )
      priorities[i] = heurdata->neighborhoods[i]->priority;

   initseed = (unsigned int)heurdata->seed + SCIPgetNVars(scip);

   /* active neighborhoods might change between init calls, reset functionality must take this into account */
   if( heurdata->bandit != NULL && SCIPbanditGetNActions(heurdata->bandit) != heurdata->nactiveneighborhoods )
   {
      SCIP_CALL( SCIPfreeBandit(scip, &heurdata->bandit) );

      heurdata->bandit = NULL;
   }

   if( heurdata->nactiveneighborhoods > 0 )
   {  /* create or reset bandit algorithm */
      if( heurdata->bandit == NULL )
      {
         SCIP_CALL( createBandit(scip, heurdata, priorities, initseed) );

         resetMinimumImprovement(heurdata);
         resetTargetNodeLimit(heurdata);
      }
      else if( heurdata->resetweights )
      {
         SCIP_CALL( SCIPresetBandit(scip, heurdata->bandit, priorities, initseed) );

         resetMinimumImprovement(heurdata);
         resetTargetNodeLimit(heurdata);
      }
   }

   heurdata->usednodes = 0;
   heurdata->ninitneighborhoods = heurdata->nactiveneighborhoods;
   resetCurrentNeighborhood(heurdata);

   SCIPfreeBufferArray(scip, &priorities);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitAlns)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free neighborhood specific data */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      NH* neighborhood = heurdata->neighborhoods[i];

      SCIP_CALL( neighborhoodExit(scip, neighborhood) );
   }

   if( heurdata->rewardfile != NULL )
   {
      fclose(heurdata->rewardfile);
      heurdata->rewardfile = NULL;
   }

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeAlns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* bandits are only initialized if a problem has been read */
   if( heurdata->bandit != NULL )
   {
      SCIP_CALL( SCIPfreeBandit(scip, &heurdata->bandit) );
   }

   /* free neighborhoods */
   for( i = 0; i < heurdata->nneighborhoods; ++i )
   {
      SCIP_CALL( alnsFreeNeighborhood(scip, &(heurdata->neighborhoods[i])) );
   }

   SCIPfreeBlockMemoryArray(scip, &heurdata->neighborhoods, NNEIGHBORHOODS);

   SCIPfreeBlockMemory(scip, &heurdata);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputNeighborhood)
{ /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(SCIPfindHeur(scip, HEUR_NAME) != NULL);
   heurdata = SCIPheurGetData(SCIPfindHeur(scip, HEUR_NAME));
   assert(heurdata != NULL);

   printNeighborhoodStatistics(scip, heurdata, file);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the alns primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurAlns(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create alns primal heuristic data */
   heurdata = NULL;
   heur = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* TODO make this a user parameter? */
   heurdata->lplimfac = LPLIMFAC;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->neighborhoods, NNEIGHBORHOODS) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecAlns, heurdata) );

   assert(heur != NULL);

   /* include all neighborhoods */
   SCIP_CALL( includeNeighborhoods(scip, heurdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyAlns) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeAlns) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitAlns) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolAlns) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitAlns) );

   /* add alns primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "offset added to the nodes budget",
         &heurdata->nodesoffset, FALSE, DEFAULT_NODESOFFSET, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start a sub-SCIP",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/waitingnodes",
         "number of nodes since last incumbent solution that the heuristic should wait",
         &heurdata->waitingnodes, TRUE, DEFAULT_WAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "fraction of nodes compared to the main SCIP for budget computation",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/startminimprove",
         "initial factor by which ALNS should at least improve the incumbent",
         &heurdata->startminimprove, TRUE, DEFAULT_STARTMINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprovelow",
         "lower threshold for the minimal improvement over the incumbent",
         &heurdata->minimprovelow, TRUE, DEFAULT_MINIMPROVELOW, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprovehigh",
         "upper bound for the minimal improvement over the incumbent",
         &heurdata->minimprovehigh, TRUE, DEFAULT_MINIMPROVEHIGH, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nsolslim",
         "limit on the number of improving solutions in a sub-SCIP call",
         &heurdata->nsolslim, FALSE, DEFAULT_NSOLSLIM, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/banditalgo",
         "the bandit algorithm: (u)pper confidence bounds, (e)xp.3, epsilon (g)reedy",
         &heurdata->banditalgo, TRUE, DEFAULT_BANDITALGO, "ueg", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/gamma",
         "weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution for exp3",
         &heurdata->exp3_gamma, TRUE, DEFAULT_GAMMA, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/beta",
         "reward offset between 0 and 1 at every observation for Exp.3",
         &heurdata->exp3_beta, TRUE, DEFAULT_BETA, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/alpha",
            "parameter to increase the confidence width in UCB",
            &heurdata->ucb_alpha, TRUE, DEFAULT_ALPHA, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usedistances",
         "distances from fixed variables be used for variable prioritization",
         &heurdata->usedistances, TRUE, DEFAULT_USEDISTANCES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useredcost",
         "should reduced cost scores be used for variable prioritization?",
         &heurdata->useredcost, TRUE, DEFAULT_USEREDCOST, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/domorefixings",
         "should the ALNS heuristic do more fixings by itself based on variable prioritization"
         "until the target fixing rate is reached?",
         &heurdata->domorefixings, TRUE, DEFAULT_DOMOREFIXINGS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/adjustfixingrate",
         "should the heuristic adjust the target fixing rate based on the success?",
         &heurdata->adjustfixingrate, TRUE, DEFAULT_ADJUSTFIXINGRATE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usesubscipheurs",
         "should the heuristic activate other sub-SCIP heuristics during its search?",
         &heurdata->usesubscipheurs, TRUE, DEFAULT_USESUBSCIPHEURS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/rewardcontrol",
         "reward control to increase the weight of the simple solution indicator and decrease the weight of the closed gap reward",
         &heurdata->rewardcontrol, TRUE, DEFAULT_REWARDCONTROL, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/targetnodefactor",
         "factor by which target node number is increased/decreased at every adjustment",
         &heurdata->targetnodefactor, TRUE, DEFAULT_TARGETNODEFACTOR, 1.0, 1e+5, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/stallnodefactor",
         "stall node limit as a fraction of total node limit",
         &heurdata->stallnodefactor, TRUE, DEFAULT_STALLNODEFACTOR, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/seed",
         "initial random seed for bandit algorithms and random decisions by neighborhoods",
         &heurdata->seed, FALSE, DEFAULT_SEED, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/adjustminimprove",
         "should the factor by which the minimum improvement is bound be dynamically updated?",
         &heurdata->adjustminimprove, TRUE, DEFAULT_ADJUSTMINIMPROVE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/eps",
         "increase exploration in epsilon-greedy bandit algorithm",
         &heurdata->epsgreedy_eps, TRUE, DEFAULT_EPS, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/rewardbaseline",
         "the reward baseline to separate successful and failed calls",
         &heurdata->rewardbaseline, TRUE, DEFAULT_REWARDBASELINE, 0.0, 0.99, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/resetweights",
         "should the bandit algorithms be reset when a new problem is read?",
         &heurdata->resetweights, TRUE, DEFAULT_RESETWEIGHTS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/" HEUR_NAME "/rewardfilename", "file name to store all rewards and the selection of the bandit",
         &heurdata->rewardfilename, TRUE, DEFAULT_REWARDFILENAME, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/subsciprandseeds",
         "should random seeds of sub-SCIPs be altered to increase diversification?",
         &heurdata->subsciprandseeds, TRUE, DEFAULT_SUBSCIPRANDSEEDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/scalebyeffort",
         "should the reward be scaled by the effort?",
         &heurdata->scalebyeffort, TRUE, DEFAULT_SCALEBYEFFORT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should cutting planes be copied to the sub-SCIP?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/fixtol",
         "tolerance by which the fixing rate may be missed/exceeded without generic (unfixing)",
         &heurdata->fixtol, TRUE, DEFAULT_FIXTOL, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselocalredcost",
         "should local reduced costs be used for generic (un)fixing?",
         &heurdata->uselocalredcost, TRUE, DEFAULT_USELOCALREDCOST, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usepscost",
            "should pseudo cost scores be used for variable priorization?",
            &heurdata->usepscost, TRUE, DEFAULT_USEPSCOST, NULL, NULL) );

   assert(SCIPfindTable(scip, TABLE_NAME_NEIGHBORHOOD) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_NEIGHBORHOOD, TABLE_DESC_NEIGHBORHOOD, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputNeighborhood,
         NULL, TABLE_POSITION_NEIGHBORHOOD, TABLE_EARLIEST_STAGE_NEIGHBORHOOD) );

   return SCIP_OKAY;
}
