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

/**@file   event_estim.c
 * @brief  event handler for tree size estimation and restarts
 *
 * This event handler plugin provides different methods for approximating the current fraction of the search
 * that has already been completed and for estimating the total tree size at completion.
 * It can trigger restarts of the current run if the current run seems hopeless.
 *
 * For details about the available approximations of search completion, please see
 *
 * Anderson, Hendel, Le Bodic, Pfetsch
 * Estimating The Size of Branch-and-Bound Trees
 * under preparation
 *
 * This code is a largely enriched version of a code that was used for clairvoyant restarts, see
 *
 * Anderson, Hendel, Le Bodic, Viernickel
 * Clairvoyant Restarts in Branch-and-Bound Search Using Online Tree-Size Estimation
 * AAAI-19: Proceedings of the Thirty-Third AAAI Conference on Artificial Intelligence, 2018
 *
 * @author Gregor Hendel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "blockmemshell/memory.h"
#include "scip/event_estim.h"
#include "scip/prop_symmetry.h"
#include "scip/pub_disp.h"
#include "scip/pub_event.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_tree.h"
#include "scip/scip_disp.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_pricer.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_table.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/type_disp.h"
#include "scip/type_event.h"
#include "scip/type_message.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_stat.h"
#include "scip/type_table.h"

#define EVENTHDLR_NAME         "estim"
#define EVENTHDLR_DESC         "event handler for tree size estimation and restarts"
#define EVENTTYPE_ESTIM      (SCIP_EVENTTYPE_NODEDELETE | SCIP_EVENTTYPE_NODEBRANCHED)

/*
 * Data structures
 */

/** enumerator for available restart policies */
enum RestartPolicy
{
   RESTARTPOLICY_NEVER      = 0,             /**< never restart (disable this event handler) */
   RESTARTPOLICY_ALWAYS     = 1,             /**< always restart (can be fine tuned by using minimum number of nodes and restart limit) */
   RESTARTPOLICY_ESTIMATION = 2,             /**< base restart on the estimation method */
   RESTARTPOLICY_COMPLETION = 3              /**< trigger restart based on search completion approximation */
};

typedef enum RestartPolicy RESTARTPOLICY;

#define RESTARTPOLICY_CHAR_NEVER        'n'
#define RESTARTPOLICY_CHAR_ALWAYS       'a'
#define RESTARTPOLICY_CHAR_COMPLETION   'c'
#define RESTARTPOLICY_CHAR_ESTIMATION   'e'

#define DES_USETRENDINLEVEL              TRUE /**< Should the trend be used in the level update? */

/* constants for the table estimation */
#define TABLE_NAME              "estim"
#define TABLE_DESC              "tree size estimations statistics table"
#define TABLE_POSITION          18500           /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE    SCIP_STAGE_INIT /**< output of the statistics table is only printed from this stage onwards */

/* constants for the search completion display column */
#define DISP_NAME               "completed"
#define DISP_DESC               "completion of search in percent (based on tree size estimation)"
#define DISP_HEADER             "compl."
#define DISP_WIDTH              8       /**< the width of the display column */
#define DISP_PRIORITY           110000  /**< the priority of the display column */
#define DISP_POSITION           30100   /**< the relative position of the display column */
#define DISP_STRIPLINE          TRUE    /**< the default for whether the display column should be separated
                                         *   with a line from its right neighbor */
#define INITIALSIZE             100
#define SESCOEFF                0.75    /**< coefficient of single exponential smoothing of estimation */

/* double exponential smoothing parameters for different time series */
#define DES_ALPHA_TREEWEIGHT 0.65
#define DES_BETA_TREEWEIGHT  0.15

#define DES_ALPHA_GAP 0.6
#define DES_BETA_GAP  0.15

#define DES_ALPHA_LEAFFREQUENCY 0.3
#define DES_BETA_LEAFFREQUENCY  0.33

#define DES_ALPHA_SSG 0.6
#define DES_BETA_SSG  0.15

#define DES_ALPHA_OPENNODES 0.6
#define DES_BETA_OPENNODES  0.15

#define MAX_REGFORESTSIZE 10000000               /**< size limit (number of nodes) for regression forest */



/* computation of search completion */
#define COMPLETIONTYPE_AUTO      'a'             /**< automatic (regression forest if available, else monotone regression on binary and SSG on nonbinary trees) */
#define COMPLETIONTYPE_REGFOREST 'r'             /**< regression forest (must be provided by user) */
#define COMPLETIONTYPE_MONOREG   'm'             /**< monotone regression (using tree weight and SSG) */
#define COMPLETIONTYPE_TREEWEIGHT 'w'            /**< use tree weight value as approximation of search tree completion */
#define COMPLETIONTYPE_SSG       's'             /**< use SSG value as approximation of search tree completion */
#define COMPLETIONTYPE_GAP       'g'             /**< use gap value as approximation of search tree completion */


/* tree size estimation method */
#define ESTIMMETHOD_COMPL        'c'             /**< estimation based on projection of current search completion */
#define ESTIMMETHOD_WBE          'b'             /**< weighted backtrack estimation */
#define ESTIMMETHOD_ENSMBL       'e'             /**< estimation based on an ensemble of the individual estimations */
#define ESTIMMETHOD_GAP          'g'             /**< estimation based on double exponential smoothing for open nodes */
#define ESTIMMETHOD_LFREQ        'l'             /**< estimation based on double exponential smoothing for leaf frequency */
#define ESTIMMETHOD_OPEN         'o'             /**< estimation based on double exponential smoothing for open nodes */
#define ESTIMMETHOD_SSG          's'             /**< estimation based on double exponential smoothing for sum of subtree gaps */
#define ESTIMMETHOD_TPROF        't'             /**< estimation based on tree profile method */
#define ESTIMMETHOD_TREEWEIGHT   'w'             /**< estimation based on double exponential smoothing for tree weight */

#define ESTIMMETHODS "bceglostw"

/* constants and default values for treeprofile parameters */
#define TREEPROFILE_MINSIZE    512                /**< minimum size (depth) that tree profile can hold */
#define SSG_STARTPRIMBOUND  SCIP_INVALID          /**< initial value of primal bound used within SSG */

/** double exponential smoothing data structure */
struct DoubleExpSmooth
{
   SCIP_Real             alpha;              /**< level smoothing constant */
   SCIP_Real             beta;               /**< trend smoothing constant */
   SCIP_Real             level;              /**< estimation of the current level used for smoothing */
   SCIP_Real             trend;              /**< estimation of the current trend (slope) */
   SCIP_Real             initialvalue;       /**< the level value at 0 observations */
   SCIP_Bool             usetrendinlevel;    /**< Should the trend be used in the level update? */
   int                   n;                  /**< number of observations */
};
typedef struct DoubleExpSmooth DOUBLEEXPSMOOTH;

/** time series data structure for leaf time series
 *
 *  These time series are the basic ingredient for tree size estimation via forecasting.
 *
 *  This general class represents concrete time series such as the closed gap, tree weight, and leaf frequency.
 *  Through callbacks for data (de-)initialization and value queries, it provides a common interface
 *  to which double exponential smoothing or window forecasts can be applied.
 */
typedef struct TimeSeries TIMESERIES;

/** data structure for convenient access of tree information */
typedef struct TreeData TREEDATA;


#define NTIMESERIES 5

/** time series position in event handler time series array */
enum TsPos
{
   TSPOS_NONE         = -1,                  /**< invalid array position */
   TSPOS_GAP          =  0,                  /**< time series position of gap */
   TSPOS_TREEWEIGHT   =  1,                  /**< time series position of tree weight */
   TSPOS_LFREQ        =  2,                  /**< time series position of leaf frequency */
   TSPOS_SSG          =  3,                  /**< time series position of SSG */
   TSPOS_OPEN         =  4                   /**< time series position of open nodes */
};

typedef enum TsPos TSPOS;

/** regression forest data structure */
typedef struct SCIP_RegForest SCIP_REGFOREST;

/** statistics collected from profile used for prediction */
struct TreeProfileStats
{
   int                   maxdepth;           /**< maximum node depth encountered */
   int                   lastfulldepth;      /**< deepest layer for which all nodes have been explored */
   int                   minwaistdepth;      /**< minimum depth of the waist, i.e. the widest part of the tree */
   int                   maxwaistdepth;      /**< maximum depth of the waist, i.e. the widest part of the tree */
};

typedef struct TreeProfileStats TREEPROFILESTATS;


/** profile data structure for tree */
struct TreeProfile
{
   SCIP_Longint*         profile;            /**< array to store the tree profile */
   int                   profilesize;        /**< size of the profile array */
   TREEPROFILESTATS      stats;              /**< statistics collected from profile used for prediction */
   SCIP_Real             lastestimate;       /**< the last estimate predicted by predictTotalSizeTreeprofile() */
   TREEPROFILESTATS      lastestimatestats;  /**< tree profile statistics at last estimation */
};

typedef struct TreeProfile TREEPROFILE;

/* default values of user parameters */
#define DEFAULT_USELEAFTS            TRUE    /**< Use leaf nodes as basic observations for time series, or all nodes? */
#define DEFAULT_REPORTFREQ           -1      /**< report frequency on estimation: -1: never, 0: always, k >= 1: k times evenly during search */
#define DEFAULT_REGFORESTFILENAME    "-"     /**< default file name of user regression forest in RFCSV format */
#define DEFAULT_COEFMONOWEIGHT       0.3667  /**< coefficient of tree weight in monotone approximation of search completion */
#define DEFAULT_COEFMONOSSG          0.6333  /**< coefficient of 1 - SSG in monotone approximation of search completion */
#define DEFAULT_COMPLETIONTYPE       COMPLETIONTYPE_AUTO /**< default computation of search tree completion */
#define DEFAULT_ESTIMMETHOD          ESTIMMETHOD_TREEWEIGHT    /**< default tree size estimation method: (c)ompletion, (e)nsemble, time series forecasts on either
                                                            * (g)ap, (l)eaf frequency, (o)open nodes,
                                                            * tree (w)eight, (s)sg, or (t)ree profile or w(b)e */
#define DEFAULT_TREEPROFILE_ENABLED  FALSE   /**< Should the event handler collect data? */
#define DEFAULT_TREEPROFILE_MINNODESPERDEPTH 20.0 /**< minimum average number of nodes at each depth before producing estimations */
#define DEFAULT_RESTARTPOLICY        'e'     /**< default restart policy: (a)lways, (c)ompletion, (e)stimation, (n)ever */
#define DEFAULT_RESTARTLIMIT          1      /**< default restart limit */
#define DEFAULT_MINNODES             1000L   /**< minimum number of nodes before restart */
#define DEFAULT_COUNTONLYLEAVES      FALSE   /**< should only leaves count for the minnodes parameter? */
#define DEFAULT_RESTARTFACTOR        50.0    /**< factor by which the estimated number of nodes should exceed the current number of nodes */
#define DEFAULT_RESTARTNONLINEAR     FALSE   /**< whether to apply a restart when nonlinear constraints are present */
#define DEFAULT_RESTARTACTPRICERS    FALSE   /**< whether to apply a restart when active pricers are used */
#define DEFAULT_HITCOUNTERLIM        50      /**< limit on the number of successive samples to really trigger a restart */
#define DEFAULT_SSG_NMAXSUBTREES     -1      /**< the maximum number of individual SSG subtrees; the old split is kept if
                                               *  a new split exceeds this number of subtrees ; -1: no limit */
#define DEFAULT_SSG_NMINNODESLASTSPLIT   0L  /**< minimum number of nodes to process between two consecutive SSG splits */
#define DEFAULT_SHOWSTATS            FALSE   /**< should statistics be shown at the end? */

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_REGFOREST*       regforest;          /**< regression forest data structure */
   TIMESERIES*           timeseries[NTIMESERIES]; /**< array of time series slots */
   TREEDATA*             treedata;           /**< tree data */
   TREEPROFILE*          treeprofile;        /**< tree profile data structure */
   char*                 regforestfilename;  /**< file name of user regression forest in RFCSV format */
   SCIP_Real             restartfactor;      /**< factor by which the estimated number of nodes should exceed the current number of nodes */
   SCIP_Real             weightlastreport;   /**< tree weight at which last report was printed */
   SCIP_Real             treeprofile_minnodesperdepth;/**< minimum average number of nodes at each depth before producing estimations */
   SCIP_Real             coefmonoweight;     /**< coefficient of tree weight in monotone approximation of search completion */
   SCIP_Real             coefmonossg;        /**< coefficient of 1 - SSG in monotone approximation of search completion */
   SCIP_Longint          minnodes;           /**< minimum number of nodes in a run before restart is triggered */
   int                   restartlimit;       /**< How often should a restart be triggered? (-1 for no limit) */
   int                   nrestartsperformed; /**< number of restarts performed so far */
   int                   restarthitcounter;  /**< the number of successive samples that would trigger a restart */
   int                   hitcounterlim;      /**< limit on the number of successive samples to really trigger a restart */
   int                   nreports;           /**< the number of reports already printed */
   int                   reportfreq;         /**< report frequency on estimation: -1: never, 0:always, k >= 1: k times evenly during search */
   int                   lastrestartrun;     /**< the last run at which this event handler triggered restart */
   char                  restartpolicyparam; /**< restart policy parameter */
   char                  estimmethod;        /**< tree size estimation method: (c)ompletion, (e)nsemble, time series forecasts on either
                                               * (g)ap, (l)eaf frequency, (o)open nodes,
                                               * tree (w)eight, (s)sg, or (t)ree profile or w(b)e */
   char                  completiontypeparam;/**< approximation of search tree completion:
                                              *   (a)uto, (g)ap, tree (w)eight, (m)onotone regression, (r)egression forest, (s)sg */
   SCIP_Bool             countonlyleaves;    /**< Should only leaves count for the minnodes parameter? */
   SCIP_Bool             useleafts;          /**< Use leaf nodes as basic observations for time series, or all nodes? */
   SCIP_Bool             treeprofile_enabled;/**< Should the event handler collect treeprofile data? */
   SCIP_Bool             treeisbinary;       /**< internal flag if all branching decisions produced 2 children */
   SCIP_Bool             restartnonlinear;   /**< whether to apply a restart when nonlinear constraints are present */
   SCIP_Bool             restartactpricers;  /**< whether to apply a restart when active pricers are used */
   SCIP_Bool             showstats;          /**< should statistics be shown at the end? */
};

typedef struct SubtreeSumGap SUBTREESUMGAP;

struct TreeData
{
   SCIP_Longint          nnodes;             /**< the total number of nodes */
   SCIP_Longint          nopen;              /**< the current number of open nodes */
   SCIP_Longint          ninner;             /**< the number of inner nodes */
   SCIP_Longint          nleaves;            /**< the number of final leaf nodes */
   SCIP_Longint          nvisited;           /**< the number of visited nodes */
   long double           weight;             /**< the current tree weight (sum of leaf weights) */
   SUBTREESUMGAP*        ssg;                /**< subtree sum gap data structure */
};

struct SubtreeSumGap
{
   SCIP_Real             value;              /**< the current subtree sum gap */
   SCIP_HASHMAP*         nodes2info;         /**< map between nodes and their subtree indices */
   SCIP_PQUEUE**         subtreepqueues;     /**< array of priority queues, one for each subtree */
   SCIP_Real             scalingfactor;      /**< the current scaling factor */
   SCIP_Real             pblastsplit;        /**< primal bound when last split occurred */
   SCIP_Longint          nodelastsplit;      /**< last node at which a subtree split occurred */
   SCIP_Longint          nminnodeslastsplit; /**< minimum number of nodes to process between two consecutive SSG splits */
   int                   nmaxsubtrees;       /**< the maximum number of individual SSG subtrees; the old split is kept if
                                               *  a new split exceeds this number of subtrees ; -1: no limit */
   int                   nsubtrees;          /**< the current number n of subtrees labeled 0 .. n - 1 */
};

/** update callback of time series */
#define DECL_TIMESERIESUPDATE(x) SCIP_RETCODE x (\
   SCIP*                 scip,                   \
   TIMESERIES*           ts,                     \
   TREEDATA*             treedata,               \
   SCIP_Real*            value                   \
   )

/** time series data structure for leaf time series */
struct TimeSeries
{
   DOUBLEEXPSMOOTH       des;                /**< double exponential smoothing data structure */
   char*                 name;               /**< name of this time series */
   SCIP_Real*            vals;               /**< value array of this time series */
   SCIP_Real*            estimation;         /**< array of estimations of this time series */
   SCIP_Real             smoothestimation;   /**< smoothened estimation value */
   SCIP_Real             targetvalue;        /**< target value of this time series */
   SCIP_Real             currentvalue;       /**< current value of time series */
   SCIP_Real             initialvalue;       /**< the initial value of time series */
   SCIP_Longint          nobs;               /**< total number of observations */
   int                   valssize;           /**< size of value array */
   int                   nvals;              /**< number of values */
   int                   resolution;         /**< current (inverse of) resolution */
   SCIP_Bool             useleafts;          /**< Should this time series be recorded at leaf nodes, or at every node? */
   DECL_TIMESERIESUPDATE((*timeseriesupdate));/**< update callback at nodes */
};

/** extended node information for SSG priority queue */
struct NodeInfo
{
   SCIP_NODE*            node;               /**< search tree node */
   SCIP_Real             lowerbound;         /**< lower bound of the node at insertion into priority queue */
   int                   pos;                /**< position of this node in priority queue */
   int                   subtreeidx;         /**< subtree index of this node */
};
typedef struct NodeInfo NODEINFO;

struct SCIP_RegForest
{
   int                   ntrees;             /**< number of trees in this forest */
   int                   dim;                /**< feature dimension */
   int*                  nbegin;             /**< array of root node indices of each tree */
   int*                  child;              /**< child index pair of each internal node, or (-1, -1) for leaves */
   int*                  splitidx;           /**< data index for split at node, or -1 at a leaf */
   SCIP_Real*            value;              /**< split position at internal nodes, prediction at leaves */
   int                   size;               /**< length of node arrays */
};

/*
 * Local methods
 */

/** convert number to string and treat SCIP_INVALID as '-' */
static
char* real2String(
   SCIP_Real             num,                /**< number to convert to string */
   char*                 buf,                /**< string buffer */
   int                   digits              /**< number of decimal digits */
   )
{
   if( num == SCIP_INVALID )/*lint !e777*/
      (void) SCIPsnprintf(buf, 1, "-");
   else if( num >= 1e+20 ) /*lint !e777*/
      (void) SCIPsnprintf(buf, 3, "inf");
   else
      (void) SCIPsnprintf(buf, SCIP_MAXSTRLEN, "%10.*f", digits, num);

   return buf;
}

/** free a regression forest data structure */
static
void SCIPregForestFree(
   SCIP_REGFOREST**      regforest           /**< regression forest data structure */
   )
{
   SCIP_REGFOREST* regforestptr;

   assert(regforest != NULL);

   if( *regforest == NULL )
      return;
   regforestptr = *regforest;

   BMSfreeMemoryArrayNull(&regforestptr->nbegin);
   BMSfreeMemoryArrayNull(&regforestptr->child);
   BMSfreeMemoryArrayNull(&regforestptr->splitidx);
   BMSfreeMemoryArrayNull(&regforestptr->value);

   BMSfreeMemory(regforest);
}

/** make a prediction with regression forest */
static
SCIP_Real SCIPregForestPredict(
   SCIP_REGFOREST*       regforest,          /**< regression forest data structure */
   SCIP_Real*            datapoint           /**< a data point that matches the dimension of this regression forest */
   )
{
   int treeidx;
   SCIP_Real value = 0.0;

   assert(regforest != NULL);
   assert(datapoint != NULL);

   SCIPdebugMessage("Start prediction method of regression forest\n");

   /* loop through the trees */
   for( treeidx = 0; treeidx < regforest->ntrees; ++treeidx )
   {
      int treepos = regforest->nbegin[treeidx];
      int* childtree = &(regforest->child[2 * treepos]);
      int* splitidxtree = &(regforest->splitidx[treepos]);
      int pos = 0;
      SCIP_Real* valuetree = &(regforest->value[treepos]);

      SCIPdebugMessage("Tree %d at position %d\n", treeidx, treepos);

      /* find the correct leaf */
      while( splitidxtree[pos] != - 1 )
      {
         int goright;

         assert(splitidxtree[pos] < regforest->dim);

         goright = (datapoint[splitidxtree[pos]] > valuetree[pos]) ? 1 : 0;
         pos = childtree[2 * pos + goright];
      }

      value += valuetree[pos];
   }

   /* return the average value that the trees predict */
   return value / (SCIP_Real)(regforest->ntrees);
}

/** read a regression forest from an rfcsv file
 *
 *  TODO improve this parser to better capture wrong user input, e.g., if the dimension is wrong
 */
static
SCIP_RETCODE SCIPregForestFromFile(
   SCIP_REGFOREST**      regforest,          /**< regression forest data structure */
   const char*           filename            /**< name of file with the regression forest data */
   )
{
   SCIP_RETCODE retcode = SCIP_OKAY;
   SCIP_FILE* file;
   SCIP_REGFOREST* regforestptr;
   char buffer[SCIP_MAXSTRLEN];
   char firstlineformat[SCIP_MAXSTRLEN];
   char dataformat[SCIP_MAXSTRLEN];
   char valuestr[SCIP_MAXSTRLEN];
   SCIP_Bool error = FALSE;
   int ntrees;
   int dim;
   int size;
   int sscanret;
   int pos;
   int treepos;

   /* try to open file */
   file = SCIPfopen(filename, "r");

   if( file == NULL )
      return SCIP_NOFILE;

   /* parse read the first line that contains the number of trees, feature dimension, and total number of nodes */
   (void) SCIPsnprintf(firstlineformat, SCIP_MAXSTRLEN, "### NTREES=%%10d FEATURE_DIM=%%10d LENGTH=%%10d\n");
   if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
   {
      error = TRUE;
      SCIPerrorMessage("Could not read first line of regression file '%s'\n", filename);
      goto CLOSEFILE;
   }

   /* coverity[secure_coding] */
   sscanret = sscanf(buffer, firstlineformat, &ntrees, &dim, &size);

   if( sscanret != 3 )
   {
      error = TRUE;
      SCIPerrorMessage("Could not extract tree information from buffer line [%s]\n", buffer);
      goto CLOSEFILE;
   }

   SCIPdebugMessage("Read ntrees=%d, dim=%d, size=%d (return value %d)\n", ntrees, dim, size, sscanret);

   /* check if the tree is too big, or numbers are negative */
   if( size > MAX_REGFORESTSIZE )
   {
      error = TRUE;
      SCIPerrorMessage("Requested size %d exceeds size limit %d for regression trees", size, MAX_REGFORESTSIZE);
      goto CLOSEFILE;
   }

   if( dim <= 0 || ntrees <= 0 || size <= 0 )
   {
      error = TRUE;
      SCIPerrorMessage("Cannot create regression tree with negative size, dimension, or number of trees\n");
      goto CLOSEFILE;
   }

   /* allocate memory in regression forest data structure */
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemory(regforest), FREEFOREST );
   BMSclearMemory(*regforest);
   regforestptr = *regforest;

   /* coverity[tainted_data] */
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&regforestptr->nbegin, ntrees), FREEFOREST );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&regforestptr->child, 2 * size), FREEFOREST ); /*lint !e647*/
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&regforestptr->splitidx, size), FREEFOREST );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&regforestptr->value, size), FREEFOREST );

   regforestptr->dim = dim;
   regforestptr->size = size;
   regforestptr->ntrees = ntrees;

   SCIPdebugMessage("Random Forest allocated\n");

   /* loop through the rest of the file, which contains the comma separated node data */
   (void) SCIPsnprintf(dataformat, SCIP_MAXSTRLEN, "%%10d,%%10d,%%10d,%%10d,%%%ds\n", SCIP_MAXSTRLEN);

   pos = 0;
   treepos = 0;
   while( !SCIPfeof(file) && !error )
   {
      int node;
      char* endptr;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;

      sscanret = sscanf(buffer, dataformat,
         &node,
         &regforestptr->child[2 * pos],
         &regforestptr->child[2 * pos + 1],
         &regforestptr->splitidx[pos],
         valuestr);

      if( sscanret != 5 )
      {
         SCIPerrorMessage("Something wrong with line %d '%s'", pos + 1, buffer);
         error = TRUE;
      }

      (void)SCIPstrToRealValue(valuestr, &regforestptr->value[pos], &endptr);

      /* new root node - increase the tree index position */
      if( node == 0 )
      {
         assert(treepos < regforestptr->ntrees);

         regforestptr->nbegin[treepos++] = pos;
      }

      ++pos;
   }

   goto CLOSEFILE;

/* insufficient memory for allocating regression forest */
FREEFOREST:
   assert(retcode == SCIP_NOMEMORY);
   SCIPregForestFree(regforest);

CLOSEFILE:
   SCIPfclose(file);

   if( error )
      retcode = SCIP_INVALIDDATA;

   return retcode;
}

/** compare two tree profile statistics for equality */
static
SCIP_Bool isEqualTreeProfileStats(
   TREEPROFILESTATS*     stats,              /**< first tree profile statistics */
   TREEPROFILESTATS*     other               /**< other tree profile statistics */
   )
{
   assert(stats != NULL);
   assert(other != NULL);

   return  stats->maxdepth == other->maxdepth &&
      stats->lastfulldepth == other->lastfulldepth &&
      stats->minwaistdepth == other->minwaistdepth &&
      stats->maxwaistdepth == other->maxwaistdepth;
}

/** copy source tree profile into destination */
static
void copyTreeProfileStats(
   TREEPROFILESTATS*     dest,               /**< destination tree profile statistics */
   TREEPROFILESTATS*     src                 /**< source tree profile statistics */
   )
{
   assert(dest != NULL);
   assert(src != NULL);

   dest->maxdepth = src->maxdepth;
   dest->lastfulldepth = src->lastfulldepth;
   dest->minwaistdepth = src->minwaistdepth;
   dest->maxwaistdepth = src->maxwaistdepth;
}

/** reset tree profile statistics */
static
void resetTreeProfileStats(
   TREEPROFILESTATS*     treeprofilestats    /**< tree profile statistics */
   )
{
   assert(treeprofilestats != NULL);

   BMSclearMemory(treeprofilestats);
}


/** extend tree profile to deeper tree */
static
SCIP_RETCODE extendMemoryTreeProfile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE*          treeprofile,        /**< tree profile data structure */
   int                   mindepth            /**< minimum depth that the tree profile should hold */
   )
{
   if( mindepth < treeprofile->profilesize )
      return SCIP_OKAY;

   if( treeprofile->profile == NULL )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &treeprofile->profile, mindepth) );
      treeprofile->profilesize = mindepth;
   }
   else
   {
      int newsize;
      int nnewelems;
      SCIP_Longint* newprofile;

      newsize = SCIPcalcMemGrowSize(scip, mindepth + 1);
      nnewelems = newsize - treeprofile->profilesize;
      assert(newsize > treeprofile->profilesize);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &treeprofile->profile, newsize) );
      newprofile = &treeprofile->profile[treeprofile->profilesize];
      BMSclearMemoryArray(newprofile, nnewelems);
      treeprofile->profilesize = newsize;
   }

   return SCIP_OKAY;
}

/** create a tree profile */
static
SCIP_RETCODE createTreeProfile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE**         treeprofile         /**< pointer to store tree profile data structure */
   )
{
   assert(scip != NULL);
   assert(treeprofile != NULL);

   SCIP_CALL( SCIPallocMemory(scip, treeprofile) );

   (*treeprofile)->profile = NULL;
   (*treeprofile)->profilesize = 0;
   SCIP_CALL( extendMemoryTreeProfile(scip, *treeprofile, TREEPROFILE_MINSIZE) );

   resetTreeProfileStats(&(*treeprofile)->stats);
   resetTreeProfileStats(&(*treeprofile)->lastestimatestats);

   (*treeprofile)->lastestimate = -1.0;

   return SCIP_OKAY;
}

/** free a tree profile */
static
void freeTreeProfile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE**         treeprofile         /**< pointer to tree profile data structure */
   )
{
   assert(scip != NULL);
   assert(treeprofile != NULL);

   if( *treeprofile == NULL )
      return;

   SCIPfreeMemoryArray(scip, &(*treeprofile)->profile);

   SCIPfreeMemory(scip, treeprofile);

   *treeprofile = NULL;
}

/** update tree profile */
static
SCIP_RETCODE updateTreeProfile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE*          treeprofile,        /**< tree profile data structure */
   SCIP_NODE*            node                /**< node that should be added to the profile */
   )
{
   int nodedepth;
   unsigned long nbits;
   SCIP_Longint nodedepthcnt;
   SCIP_Longint maxnodes;

   assert(scip != NULL);
   assert(node != NULL);

   if( treeprofile == NULL )
      return SCIP_OKAY;

   nodedepth = SCIPnodeGetDepth(node);
   assert(nodedepth >= 0);
   maxnodes = treeprofile->profile[treeprofile->stats.minwaistdepth];
   assert(treeprofile->stats.minwaistdepth == treeprofile->stats.maxwaistdepth ||
      maxnodes == treeprofile->profile[treeprofile->stats.maxwaistdepth]);

   /* ensure that the memory can hold at least this depth */
   SCIP_CALL( extendMemoryTreeProfile(scip, treeprofile, nodedepth) );

   nodedepthcnt = ++treeprofile->profile[nodedepth];

   /* Is this level fully explored? We assume binary branching. The first condition ensures that the bit shift operation
    * of the second condition represents a feasible power of unsigned int. The largest power of 2 representable
    * by unsigned int is 2^{8*sizeof(unsigned int) - 1}. */
   nbits = 8*sizeof(unsigned int);
   /* coverity[overflow_before_widen] */
   if( (unsigned int)nodedepth < nbits && nodedepthcnt == (1U << nodedepth) )/*lint !e647*/
   {
      SCIPdebugMsg(scip, "Level %d fully explored: %" SCIP_LONGINT_FORMAT " nodes\n", nodedepth, nodedepthcnt);

      treeprofile->stats.lastfulldepth = nodedepth;
   }

   /* update maximum depth */
   if( treeprofile->stats.maxdepth < nodedepth )
   {
      treeprofile->stats.maxdepth = nodedepth;
      SCIPdebugMsg(scip, "Maximum depth increased to %d\n", treeprofile->stats.maxdepth);
   }

   /* minimum and maximum waist now coincide */
   if( nodedepthcnt > maxnodes )
   {
      treeprofile->stats.minwaistdepth = treeprofile->stats.maxwaistdepth = nodedepth;
      SCIPdebugMsg(scip, "Updating depth of tree waist: %d (%" SCIP_LONGINT_FORMAT " nodes)\n",
         treeprofile->stats.minwaistdepth, nodedepthcnt);
   }
   else if( nodedepthcnt == maxnodes )
   {
      /* enlarge the interval in which the waist lies */
      if( treeprofile->stats.minwaistdepth > nodedepth )
         treeprofile->stats.minwaistdepth = nodedepth;
      else if( treeprofile->stats.maxwaistdepth < nodedepth )
         treeprofile->stats.maxwaistdepth = nodedepth;
   }
   assert(treeprofile->stats.minwaistdepth <= treeprofile->stats.maxwaistdepth);

   return SCIP_OKAY;
}

/** make a prediction of the total tree size based on the current tree profile */
static
SCIP_Real predictTotalSizeTreeProfile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE*          treeprofile,        /**< tree profile data structure */
   SCIP_Real             minnodesperdepth    /**< minimum number of average nodes per depth to make a prediction */
   )
{
   SCIP_Real estimate;
   SCIP_Real growthfac;
   int d;
   int waist;

   /* prediction is disabled */
   if( treeprofile == NULL )
      return -1.0;

   /* two few nodes to make a prediction */
   if( minnodesperdepth * treeprofile->stats.maxdepth > SCIPgetNNodes(scip) )
      return -1.0;

   /* reuse previous estimation if tree profile hasn't changed */
   if( isEqualTreeProfileStats(&treeprofile->lastestimatestats, &treeprofile->stats) )
   {
      SCIPdebugMsg(scip, "Reusing previous estimation result %g\n", treeprofile->lastestimate);

      return treeprofile->lastestimate;
   }

   /* compute the (depth of the) waist as convex combination between the minimum and maximum waist depths */
   waist = (2 * treeprofile->stats.maxwaistdepth + treeprofile->stats.minwaistdepth) / 3;

   growthfac = 2;
   estimate = 1;

   /* loop over all full levels */
   for( d = 1; d < treeprofile->stats.lastfulldepth; ++d )
   {
      SCIP_Real gamma_d = 2.0;

      estimate += growthfac;
      growthfac *= gamma_d;
   }

   /* loop until the waist is reached */
   for( ; d < waist; ++d )
   {
      SCIP_Real gamma_d = 2.0 - (d - treeprofile->stats.lastfulldepth + 1.0)/(waist - treeprofile->stats.lastfulldepth + 1.0);

      assert(1.0 <= gamma_d && gamma_d <= 2.0);
      estimate += growthfac;
      growthfac *= gamma_d;
   }

   /* loop over the remaining levels */
   for( ; d <= treeprofile->stats.maxdepth; ++d )
   {
      SCIP_Real gamma_d = (1.0 - (d - waist + 1.0)/(treeprofile->stats.maxdepth - waist + 1.0));
      assert(0.0 <= gamma_d && gamma_d <= 1.0);

      estimate += growthfac;
      growthfac *= gamma_d;
   }

   /* copy tree profile statistics */
   copyTreeProfileStats(&treeprofile->lastestimatestats, &treeprofile->stats);

   treeprofile->lastestimate = estimate;

   return estimate;
}

/** clean subtrees stored as priority queues */
static
void subtreeSumGapDelSubtrees(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg                 /**< subtree sum gap data structure */
   )
{
   assert(ssg->nsubtrees <= 1 || ssg->subtreepqueues != NULL);

   /* free all previous priority queues */
   if( ssg->nsubtrees > 1 )
   {
      int s;

      for( s = 0; s < ssg->nsubtrees; ++s )
      {
         int i;
         SCIP_PQUEUE* pqueue = ssg->subtreepqueues[s];
         NODEINFO** nodeinfos;

         assert(pqueue != NULL);
         nodeinfos = (NODEINFO**)SCIPpqueueElems(pqueue);

         /* free all remaining elements in reverse order */
         for( i = SCIPpqueueNElems(pqueue); --i >= 0; )
         {
            NODEINFO* nodeinfo = nodeinfos[i];
            assert(nodeinfo != NULL);
            SCIPfreeBlockMemory(scip, &nodeinfo);
         }

         SCIPpqueueFree(&pqueue);
      }

      SCIPfreeBlockMemoryArray(scip, &ssg->subtreepqueues, ssg->nsubtrees);
   }

   ssg->subtreepqueues = NULL;
}

/** reset subtree sum gap */
static
SCIP_RETCODE subtreeSumGapReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg                 /**< subtree sum gap data structure */
   )
{
   assert(ssg != NULL);
   assert(ssg->nodes2info != NULL);

   SCIP_CALL( SCIPhashmapRemoveAll(ssg->nodes2info) );

   subtreeSumGapDelSubtrees(scip, ssg);

   ssg->value = 1.0;
   ssg->scalingfactor = 1.0;
   ssg->nsubtrees = 1;
   ssg->subtreepqueues = NULL;
   ssg->pblastsplit = SSG_STARTPRIMBOUND;
   ssg->nodelastsplit = -1L;

   return SCIP_OKAY;
}

/** create a subtree sum gap */
static
SCIP_RETCODE subtreeSumGapCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP**       ssg                 /**< pointer to store subtree sum gap data structure */
   )
{
   assert(scip != NULL);
   assert(ssg != NULL);

   /* allocate storage */
   SCIP_CALL( SCIPallocMemory(scip, ssg) );
   SCIP_CALL( SCIPhashmapCreate(&(*ssg)->nodes2info, SCIPblkmem(scip), INITIALSIZE) );

   /* explicitly set this to skip removal of subtrees during reset */
   (*ssg)->nsubtrees = 0;

   /* reset ssg */
   SCIP_CALL( subtreeSumGapReset(scip, *ssg) );

   return SCIP_OKAY;
}

/** free a subtree sum gap */
static
void subtreeSumGapFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP**       ssg                 /**< pointer to store subtree sum gap data structure */
   )
{
   assert(scip != NULL);

   if( *ssg == NULL )
      return;

   if( (*ssg)->nodes2info != NULL )
   {
      SCIPhashmapFree(&(*ssg)->nodes2info);
   }

   /* delete all subtree data */
   subtreeSumGapDelSubtrees(scip, *ssg);

   SCIPfreeMemory(scip, ssg);
}

/** compare two node infos by comparing their lower bound */
static
SCIP_DECL_SORTPTRCOMP(compareNodeInfos)
{
   NODEINFO* nodeinfo1 = (NODEINFO*)elem1;
   NODEINFO* nodeinfo2 = (NODEINFO*)elem2;

   if( nodeinfo1->lowerbound < nodeinfo2->lowerbound )
      return -1;
   else if( nodeinfo1->lowerbound > nodeinfo2->lowerbound )
      return 1;

   return 0;
}

/** position change callback of element in priority queue */
static
SCIP_DECL_PQUEUEELEMCHGPOS(elemChgPosNodeInfo)
{
   NODEINFO* nodeinfo = (NODEINFO*)elem;

   assert(oldpos == -1 || oldpos == nodeinfo->pos);
   nodeinfo->pos = newpos;
}

/** store node in SSG data structure */
static
SCIP_RETCODE subtreeSumGapStoreNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg,                /**< subtree sum gap data structure */
   SCIP_NODE*            node,               /**< node that should be stored */
   int                   subtreeidx          /**< subtree index of that node */
   )
{
   NODEINFO* nodeinfo;

   assert(scip != NULL);
   assert(ssg != NULL);
   assert(node != NULL);

   /* create a new node info */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nodeinfo) );

   /* store node information in data structure and insert into priority queue */
   nodeinfo->node = node;
   nodeinfo->subtreeidx = subtreeidx;
   nodeinfo->pos = -1;
   nodeinfo->lowerbound = SCIPnodeGetLowerbound(node);

   SCIPdebugMsg(scip, "Inserting label %d for node number %" SCIP_LONGINT_FORMAT " (%p)\n",
      subtreeidx, SCIPnodeGetNumber(node), (void*)node);

   assert(!SCIPhashmapExists(ssg->nodes2info, (void*)node));
   /* store node information in Hash Map */
   SCIP_CALL( SCIPhashmapInsert(ssg->nodes2info, (void*)node, (void*)nodeinfo) );

   /* create the corresponding priority queue, if it does not exist yet */
   assert(subtreeidx >= 0);
   assert(subtreeidx < ssg->nsubtrees);

   if( ssg->subtreepqueues[subtreeidx] == NULL )
   {
      SCIP_CALL( SCIPpqueueCreate(&ssg->subtreepqueues[subtreeidx], 5, 1.2, compareNodeInfos, elemChgPosNodeInfo) );
   }

   /* insert node and ensure that its position is up to date */
   SCIP_CALL( SCIPpqueueInsert(ssg->subtreepqueues[subtreeidx], (void*)nodeinfo) );
   assert(0 <= nodeinfo->pos);
   assert(SCIPpqueueNElems(ssg->subtreepqueues[subtreeidx]) > nodeinfo->pos);
   assert(SCIPpqueueElems(ssg->subtreepqueues[subtreeidx])[nodeinfo->pos] == (void*)nodeinfo);

   return SCIP_OKAY;
}

/** split the open nodes of the current tree */
static
SCIP_RETCODE subtreeSumGapSplit(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg,                /**< subtree sum gap data structure */
   SCIP_Bool             addfocusnode        /**< should the focus node be a subtree, too? */
   )
{
   SCIP_NODE** opennodes[3];
   int nopennodes[3];
   int label;
   int t;
   int nnewsubtrees;

   assert(scip != NULL);
   assert(ssg != NULL);

   /* query the open nodes of SCIP */
   SCIP_CALL( SCIPgetOpenNodesData(scip, &opennodes[0], &opennodes[1], &opennodes[2], &nopennodes[0], &nopennodes[1], &nopennodes[2]) );

   nnewsubtrees = nopennodes[0] + nopennodes[1] + nopennodes[2] + (addfocusnode ? 1 : 0);

   /* clear hash map from entries */
   SCIP_CALL( SCIPhashmapRemoveAll(ssg->nodes2info) );

   /* delete all subtrees */
   subtreeSumGapDelSubtrees(scip, ssg);

   ssg->nsubtrees = nnewsubtrees;
   SCIPdebugMsg(scip, "Splitting tree into %d subtrees\n", ssg->nsubtrees);

   /* create priority queue array */
   if( ssg->nsubtrees > 1 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &ssg->subtreepqueues, ssg->nsubtrees) );
   }
   else
   {
      ssg->subtreepqueues = NULL;

      return SCIP_OKAY;
   }

   /* loop over node types (leaves, siblings, children) */
   label = 0;
   for( t = 0; t < 3; ++t )
   {
      SCIP_NODE** nodes = opennodes[t];
      int nnodes = nopennodes[t];
      int n;

      /* label each open node as new, separate subtree */
      for( n = 0; n < nnodes; ++n )
      {
         SCIP_NODE* node = nodes[n];
         SCIP_CALL( subtreeSumGapStoreNode(scip, ssg, node, label++) );
      }
   }

   if( addfocusnode )
   {
      assert(SCIPgetFocusNode(scip) != NULL);
      SCIP_CALL( subtreeSumGapStoreNode(scip, ssg, SCIPgetFocusNode(scip), label) );
   }

   return SCIP_OKAY;
}

/** compute a gap between a lower bound and the current upper bound */
static
SCIP_Real calcGap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lowerbound          /**< lower bound value */
   )
{
   SCIP_Real db;
   SCIP_Real pb;
   SCIP_Real abspb;
   SCIP_Real absdb;
   SCIP_Real gap;

   if( SCIPisInfinity(scip, lowerbound) || lowerbound >= SCIPgetUpperbound(scip) )
      return 0.0;

   if( SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
      return 1.0;

   db = SCIPretransformObj(scip, lowerbound);
   pb = SCIPgetPrimalbound(scip);

   if( SCIPisEQ(scip, db, pb) )
      return 0.0;

   abspb = REALABS(pb);
   absdb = REALABS(db);
   gap = REALABS(pb - db)/MAX(abspb,absdb);
   gap = MIN(gap, 1.0);

   return gap;
}

/** remove node from the subtree sum gap (because it has been solved by branching or is a leaf) */
static
SCIP_RETCODE subtreeSumGapRemoveNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg,                /**< subtree sum gap data structure */
   SCIP_NODE*            node                /**< node that should be removed */
   )
{
   NODEINFO* nodeinfo;
   int subtreeidx;
   int pos;
   SCIP_PQUEUE* pqueue;

   if( ssg->nsubtrees <= 1 )
      return SCIP_OKAY;

   nodeinfo = (NODEINFO*)SCIPhashmapGetImage(ssg->nodes2info, (void*)node);

   /* it can happen that the node was not created via branching; search for the most recent ancestor in the queue */
   if( nodeinfo == NULL )
   {
      do
      {
         node = SCIPnodeGetParent(node);
      } while( node != NULL && (nodeinfo = (NODEINFO*)SCIPhashmapGetImage(ssg->nodes2info, (void*)node)) == NULL);

      /* no ancestor found */
      if( nodeinfo == NULL )
         return SCIP_OKAY;
   }

   /* get open nodes of this subtree stored as priority queue */
   subtreeidx = nodeinfo->subtreeidx;
   pqueue = ssg->subtreepqueues[subtreeidx];
   assert(pqueue != NULL);

   /* delete the element from the priority queue */
   pos = nodeinfo->pos;
   assert(pos >= 0);
   assert(pos < SCIPpqueueNElems(pqueue));
   assert(SCIPpqueueElems(pqueue)[pos] == (void *)nodeinfo);
   SCIPpqueueDelPos(pqueue, pos);

   /* update ssg if removed node was the lower bound defining node of its subtree */
   if( pos == 0 )
   {
      NODEINFO* nodeinfofirst;
      SCIP_Real oldgap;
      SCIP_Real newgap;

      oldgap = calcGap(scip, nodeinfo->lowerbound);
      nodeinfofirst = (NODEINFO*)SCIPpqueueFirst(ssg->subtreepqueues[subtreeidx]);
      assert(nodeinfofirst == NULL || subtreeidx == nodeinfofirst->subtreeidx);
      newgap = calcGap(scip, nodeinfofirst != NULL ? nodeinfofirst->lowerbound : SCIPinfinity(scip) );

      assert(SCIPisLE(scip, newgap, oldgap));

      /* the SSG value is always up-to-date because it is recomputed when the primal bound changes */
      ssg->value += ssg->scalingfactor * MIN(newgap - oldgap, 0.0);
   }

   SCIP_CALL( SCIPhashmapRemove(ssg->nodes2info, (void*)node) );

   SCIPdebugMsg(scip, "Removed node %" SCIP_LONGINT_FORMAT " from open nodes of SSG\n",
      SCIPnodeGetNumber(node));

   SCIPfreeBlockMemory(scip, &nodeinfo);

   return SCIP_OKAY;
}

/** insert children into subtree sum gap */
static
SCIP_RETCODE subtreeSumGapInsertChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg                 /**< subtree sum gap data structure */
   )
{
   int nchildren;
   SCIP_NODE** children;
   SCIP_NODE* focusnode;
   SCIP_NODE* parentnode;
   NODEINFO* parentnodeinfo;
   int parentnodelabel;
   int n;

   assert(scip != NULL);
   assert(ssg != NULL);

   if( ssg->nsubtrees == 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetChildren(scip, &children, &nchildren) );

   if( nchildren == 0 )
      return SCIP_OKAY;

   focusnode = SCIPgetFocusNode(scip);

   /* a rare case: the search has been stopped at some point, and the current focus node is the only descendant
    * of its parent node
    */
   if( !SCIPhashmapExists(ssg->nodes2info, (void*)focusnode) )
   {
      parentnode = focusnode;
      do
      {
         parentnode = SCIPnodeGetParent(parentnode);
      } while( parentnode != NULL && !SCIPhashmapExists(ssg->nodes2info, (void *)parentnode));

      assert(parentnode != NULL && SCIPhashmapExists(ssg->nodes2info, (void *)parentnode));
   }
   else
      parentnode = focusnode;

   parentnodeinfo = (NODEINFO*)SCIPhashmapGetImage(ssg->nodes2info, (void *)parentnode);
   parentnodelabel = parentnodeinfo->subtreeidx;

   /* loop over children and insert the focus node label */
   for( n = 0; n < nchildren; ++n )
   {
      assert(SCIPnodeGetParent(children[n]) == focusnode);

      SCIPdebugMsg(scip,
         "Inserting label %d for node number %" SCIP_LONGINT_FORMAT " (parent %" SCIP_LONGINT_FORMAT ")\n",
         parentnodelabel, SCIPnodeGetNumber(children[n]), SCIPnodeGetNumber(parentnode));

      SCIP_CALL( subtreeSumGapStoreNode(scip, ssg, children[n], parentnodelabel) );
   }

   /* remove focus node from hash map */
   SCIP_CALL( subtreeSumGapRemoveNode(scip, ssg, parentnode) );

   return SCIP_OKAY;
}

/* this function is inefficient because it loops over all open nodes, but can be used for debugging */
#ifdef SCIP_DISABLED_CODE
/** compute subtree sum gap from scratch (inefficiently because loop over all open nodes) */
static
SCIP_RETCODE subtreesumgapComputeFromScratch(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg,                /**< subtree sum gap data structure */
   SCIP_Bool             updatescaling       /**< should the scaling factor be updated? */
   )
{
   SCIP_Real* lowerbounds;
   SCIP_NODE** opennodes[3];
   SCIP_Real gapsum = 0;
   SCIP_Real pb;
   int nopennodes[3];
   int l;
   int t;

   /* treat trivial cases: only 1 subtree, no incumbent solution */
   if( SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
   {
      ssg->value = 1.0;

      return SCIP_OKAY;
   }

   /* simply use normal gap in trivial case */
   if( ssg->nsubtrees == 1 )
   {
      ssg->value = calcGap(scip, SCIPgetLowerbound(scip));

      return SCIP_OKAY;
   }

    /* allocate temporary memory to store lower bound for every subtree    */
   SCIP_CALL( SCIPallocBufferArray(scip, &lowerbounds, ssg->nsubtrees) );

    /* initialize lower bounds as SCIPinfinity(scip) */
   for( l = 0; l < ssg->nsubtrees; ++l )
      lowerbounds[l] = SCIPinfinity(scip);

    /* loop over children, siblings, and leaves to update subtree lower bounds */
   SCIP_CALL( SCIPgetOpenNodesData(scip, &opennodes[0], &opennodes[1], &opennodes[2], &nopennodes[0], &nopennodes[1], &nopennodes[2]) );

   /* loop over the three types leaves, siblings, leaves */
   for( t = 0; t < 3; ++t )
   {
      int n;
      /* loop over nodes of this type */
      for( n = 0; n < nopennodes[t]; ++n )
      {
         SCIP_NODE* node = opennodes[t][n];
         NODEINFO* nodeinfo;
         SCIP_Real lowerbound;
         int label;
         nodeinfo = (NODEINFO*)SCIPhashmapGetImage(ssg->nodes2info, (void *)node);
         label = nodeinfo->subtreeidx;
         lowerbound = nodeinfo->lowerbound;

         assert(label >= 0 && label < ssg->nsubtrees);
         lowerbounds[label] = MIN(lowerbounds[label], lowerbound);
      }
   }

   /* compute subtree gaps in original space; sum them up */
   pb = SCIPgetPrimalbound(scip);
   for( l = 0; l < ssg->nsubtrees; ++l )
   {
      SCIP_Real subtreedualbound;
      SCIP_Real subtreegap;
      /* skip subtrees with infinite lower bound; they are empty and contribute 0.0 to the gap sum term */
      if( SCIPisInfinity(scip, lowerbounds[l]) )
         continue;

      subtreedualbound = SCIPretransformObj(scip, lowerbounds[l]);

      if( SCIPisEQ(scip, subtreedualbound, pb) )
         continue;

      subtreegap = REALABS(pb - subtreedualbound)/MAX(REALABS(pb),REALABS(subtreedualbound));
      subtreegap = MIN(subtreegap, 1.0);

      gapsum += subtreegap;
   }

   /* update the scaling factor by using the previous SSG value divided by the current gapsum */
   if( updatescaling )
   {
      ssg->scalingfactor = ssg->value / MAX(gapsum, 1e-6);
   }

   /* update and store SSG value by considering scaling factor */
   ssg->value = ssg->scalingfactor * gapsum;

   SCIPfreeBufferArray(scip, &lowerbounds);

   return SCIP_OKAY;
}
#endif

/** compute subtree sum gap from scratch efficiently (linear effort in the number of subtrees) */
static
SCIP_RETCODE subtreeSumGapComputeFromScratchEfficiently(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg,                /**< subtree sum gap data structure */
   SCIP_Bool             updatescaling       /**< should the scaling factor be updated? */
   )
{
   SCIP_Real gapsum = 0.0;
   int l;

   /* treat trivial cases: only 1 subtree, no incumbent solution */
   if( SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
   {
      ssg->value = 1.0;

      return SCIP_OKAY;
   }

   if( ssg->nsubtrees == 1 )
   {
      ssg->value = calcGap(scip, SCIPgetLowerbound(scip));

      return SCIP_OKAY;
   }

   /* compute subtree gaps in original space; sum them up */
   for( l = 0; l < ssg->nsubtrees; ++l )
   {
      SCIP_Real subtreegap;
      NODEINFO* nodeinfo;

      assert(ssg->subtreepqueues[l] != NULL);

      nodeinfo = (NODEINFO*)SCIPpqueueFirst(ssg->subtreepqueues[l]);

      /* skip subtrees with infinite lower bound; they are empty and contribute 0.0 to the gap sum term */
      if( nodeinfo == NULL || SCIPisInfinity(scip, nodeinfo->lowerbound) )
         continue;

      subtreegap = calcGap(scip, nodeinfo->lowerbound);

      gapsum += subtreegap;
   }

   /* update the scaling factor by using the previous SSG value divided by the current gapsum */
   if( updatescaling )
   {
      ssg->scalingfactor = ssg->value / MAX(gapsum, 1e-6);
   }

   /* update and store SSG value by considering scaling factor */
   ssg->value = ssg->scalingfactor * gapsum;

   return SCIP_OKAY;
}

/** update the subtree sum gap after a node event (branching or deletion of a node) */
static
SCIP_RETCODE subtreeSumGapUpdate(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBTREESUMGAP*        ssg,                /**< subtree sum gap data structure */
   SCIP_NODE*            node,               /**< the corresponding node */
   int                   nchildren,          /**< number of children */
   SCIP_Longint          nsolvednodes        /**< number of solved nodes so far, used as a time stamp */
   )
{
   SCIP_Bool updatescaling = FALSE;
   SCIP_Bool insertchildren = (ssg->nsubtrees > 1 && nchildren > 0);

   /* if the instance is solved or a node is cutoff at the initsolve stage or we are unbounded, the ssg is 0 */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED || SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE || SCIPisInfinity(scip, -SCIPgetUpperbound(scip)) )
   {
      ssg->value = 0.0;

      return SCIP_OKAY;
   }

   /* make a new tree split if the primal bound has changed. */
   if( ! SCIPisInfinity(scip, SCIPgetUpperbound(scip)) && ! SCIPisEQ(scip, SCIPgetPrimalbound(scip), ssg->pblastsplit) )
   {
      int nnewsubtrees;
      SCIP_Bool addfocusnode;

      addfocusnode = SCIPgetFocusNode(scip) != NULL && SCIPgetNChildren(scip) == 0 && !SCIPwasNodeLastBranchParent(scip, SCIPgetFocusNode(scip));
      nnewsubtrees = SCIPgetNSiblings(scip) + SCIPgetNLeaves(scip) + SCIPgetNChildren(scip) + (addfocusnode ? 1 : 0);

      /* check if number of new subtrees does not exceed maximum number of subtrees; always split if no split happened, yet */
      if( ssg->nsubtrees <= 1 ||
         ((ssg->nmaxsubtrees == -1 || nnewsubtrees <= ssg->nmaxsubtrees) &&
            (nsolvednodes - ssg->nodelastsplit >= ssg->nminnodeslastsplit)) )
      {
         SCIP_CALL( subtreeSumGapSplit(scip, ssg, addfocusnode) );

         /* remember time stamp */
         ssg->nodelastsplit = nsolvednodes;
      }
      else
      {
         if( ssg->nmaxsubtrees != -1 && nnewsubtrees >= ssg->nmaxsubtrees )
         {
            SCIPdebugMsg(scip, "Keep split into %d subtrees because new split into %d subtrees exceeds limit %d\n",
               ssg->nsubtrees, nnewsubtrees, ssg->nmaxsubtrees);
         }
         else
         {
            SCIPdebugMsg(scip, "Keep split into %d subtrees from %" SCIP_LONGINT_FORMAT " nodes ago\n",
               ssg->nsubtrees, nsolvednodes - ssg->nodelastsplit);
         }

         /* no new split has happened; insert the new children to their SSG subtree */
         if( insertchildren )
         {
            SCIP_CALL( subtreeSumGapInsertChildren(scip, ssg) );
         }
      }

      ssg->pblastsplit = SCIPgetPrimalbound(scip);

      updatescaling = TRUE;

      /* compute the current SSG value from scratch */
      SCIP_CALL( subtreeSumGapComputeFromScratchEfficiently(scip, ssg, updatescaling) );
   }
   /* otherwise, if new children have been created, label them */
   else if( insertchildren )
   {
      SCIP_CALL( subtreeSumGapInsertChildren(scip, ssg) );
   }

   /* remove the node from the hash map if it is a leaf */
   if( nchildren == 0 )
   {
      SCIP_CALL( subtreeSumGapRemoveNode(scip, ssg, node) );
   }

   return SCIP_OKAY;
}

/** reset tree data */
static
SCIP_RETCODE resetTreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEDATA*             treedata            /**< tree data */
   )
{
   /* simply set everything to 0 */
   treedata->ninner = treedata->nleaves = treedata->nvisited = 0L;
   treedata->weight = 0.0;

   /* set up root node */
   treedata->nnodes = 1;
   treedata->nopen = 1;

   SCIP_CALL( subtreeSumGapReset(scip, treedata->ssg) );

   return SCIP_OKAY;
}

/** create tree data structure */
static
SCIP_RETCODE createTreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEDATA**            treedata            /**< pointer to store tree data */
   )
{
   assert(treedata != NULL);
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, treedata) );

   SCIP_CALL( subtreeSumGapCreate(scip, &(*treedata)->ssg) );

   SCIP_CALL( resetTreeData(scip, *treedata) );

   return SCIP_OKAY;
}

/** free tree data structure */
static
void freeTreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEDATA**            treedata            /**< pointer to tree data */
   )
{
   assert(scip != NULL);

   if( *treedata == NULL )
      return;

   subtreeSumGapFree(scip, &(*treedata)->ssg);

   SCIPfreeMemory(scip, treedata);
   *treedata = NULL;
}

/** update tree data structure after a node has been solved/is about to be deleted */
static
SCIP_RETCODE updateTreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEDATA*             treedata,           /**< tree data */
   SCIP_NODE*            node,               /**< the corresponding node */
   int                   nchildren           /**< the number of children */
   )
{
   assert(node != NULL);

   ++treedata->nvisited;
   treedata->nopen--;

   if( nchildren == 0 )
   {
      int depth = SCIPnodeGetDepth(node);
      treedata->nleaves++;
      treedata->weight += pow(0.5, (SCIP_Real)depth);
   }
   else
   {
      treedata->nnodes += nchildren;
      treedata->nopen += nchildren;
      ++treedata->ninner;
   }

   /* update the subtree sum gap */
   if( ! SCIPisInRestart(scip) )
   {
      SCIP_CALL( subtreeSumGapUpdate(scip, treedata->ssg, node, nchildren, treedata->nvisited) );
   }

   return SCIP_OKAY;
}

/** get weighted backtrack estimation from this tree data */
static
SCIP_Real treeDataGetWbe(
   TREEDATA*             treedata            /**< tree data */
   )
{
   if( treedata->weight <= 0.0 || treedata->nleaves == 0 )
      return -1.0;

   return 2.0 * treedata->nleaves / (SCIP_Real)treedata->weight - 1.0;
}

#ifdef SCIP_DEBUG
/* print method for tree data */
static
char* treeDataPrint(
   TREEDATA*             treedata,           /**< tree data */
   char*                 strbuf              /**< string buffer */
   )
{
   (void )SCIPsnprintf(strbuf, SCIP_MAXSTRLEN,
      "Tree Data: %" SCIP_LONGINT_FORMAT " nodes ("
      "%" SCIP_LONGINT_FORMAT " visited, "
      "%" SCIP_LONGINT_FORMAT " inner, "
      "%" SCIP_LONGINT_FORMAT " leaves, "
      "%" SCIP_LONGINT_FORMAT " open), "
      "weight: %.4Lf, ssg %.4f",
      treedata->nnodes,
      treedata->nvisited,
      treedata->ninner,
      treedata->nleaves,
      treedata->nopen,
      treedata->weight,
      treedata->ssg->value
      );
   return strbuf;
}
#endif

/** reset double exponential smoothing */
static
void doubleExpSmoothReset(
   DOUBLEEXPSMOOTH*      des,                /**< double exponential smoothing data structure */
   SCIP_Real             initialvalue        /**< the initial value */
   )
{
   des->n = 0;
   des->level = SCIP_INVALID;
   des->trend = SCIP_INVALID;
   des->initialvalue = initialvalue;
}

/** initialize a double exponential smoothing data structure */
static
void doubleExpSmoothInit(
   DOUBLEEXPSMOOTH*      des,                /**< double exponential smoothing data structure */
   SCIP_Real             x1                  /**< the first sample value */
   )
{
   assert(des != NULL);

   des->n = 1;
   des->level = x1;
   des->trend = x1 - des->initialvalue;

   des->usetrendinlevel = DES_USETRENDINLEVEL;

   return;
}

/** update a double exponential smoothing data structure */
static
void doubleExpSmoothUpdate(
   DOUBLEEXPSMOOTH*      des,                /**< double exponential smoothing data structure */
   SCIP_Real             xnew                /**< new sample value */
   )
{
   if( des->n == 0 )
      doubleExpSmoothInit(des, xnew);
   else
   {
      SCIP_Real newlevel;
      SCIP_Real newtrend;

      newlevel = des->alpha * xnew + (1.0 - des->alpha) * (des->level + des->usetrendinlevel ? des->trend : 0.0);
      newtrend = des->beta * (newlevel - des->level) + (1.0 - des->beta) * des->trend;

      des->level = newlevel;
      des->trend = newtrend;
   }
}

/** get the current trend (slope) computed by this double exponential smoothing */
static
SCIP_Real doubleExpSmoothGetTrend(
   DOUBLEEXPSMOOTH*      des                 /**< double exponential smoothing data structure */
   )
{
   assert(des != NULL);

   if( des->n == 0 )
      return SCIP_INVALID;

   return des->trend;
}

/** reset time series */
static
void timeSeriesReset(
   TIMESERIES*           timeseries          /**< pointer to store time series */
   )
{
   timeseries->resolution = 1;
   timeseries->nvals = 0;
   timeseries->nobs = 0L;
   timeseries->currentvalue = timeseries->initialvalue;
   timeseries->smoothestimation = SCIP_INVALID;

   doubleExpSmoothReset(&timeseries->des, timeseries->initialvalue);
}

/** create a time series object */
static
SCIP_RETCODE timeSeriesCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMESERIES**          timeseries,         /**< pointer to store time series */
   const char*           name,               /**< name of this time series */
   SCIP_Real             targetvalue,        /**< target value of this time series */
   SCIP_Real             initialvalue,       /**< the initial value of time series */
   SCIP_Real             alpha,              /**< alpha parameter (level weight) for double exponential smoothing */
   SCIP_Real             beta,               /**< beta parameter (level weight) for double exponential smoothing */
   DECL_TIMESERIESUPDATE ((*timeseriesupdate)) /**< update callback at nodes, or NULL */
   )
{
   TIMESERIES* timeseriesptr;
   assert(scip != NULL);
   assert(timeseries != NULL);
   assert(name != NULL);
   assert(alpha >= 0.0 && alpha <= 1);
   assert(beta >= 0.0 && beta <= 1);

   SCIP_CALL( SCIPallocMemory(scip, timeseries) );

   timeseriesptr = *timeseries;
   assert(timeseriesptr != NULL);

   /* copy name */
   SCIP_ALLOC( BMSduplicateMemoryArray(&timeseriesptr->name, name, strlen(name)+1) );

   /* copy callbacks */
   assert(timeseriesupdate != NULL);
   timeseriesptr->timeseriesupdate = timeseriesupdate;

   timeseriesptr->targetvalue = targetvalue;
   timeseriesptr->valssize = 1024;
   timeseriesptr->initialvalue = initialvalue;

   SCIP_CALL( SCIPallocMemoryArray(scip, &timeseriesptr->vals, timeseriesptr->valssize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &timeseriesptr->estimation, timeseriesptr->valssize) );

   timeSeriesReset(timeseriesptr);

   timeseriesptr->des.alpha = alpha;
   timeseriesptr->des.beta = beta;

   SCIPdebugMsg(scip, "Finished creation of time series '%s'\n", timeseriesptr->name);

   return SCIP_OKAY;
}

/** free a time series */
static
void timeSeriesFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMESERIES**          timeseries          /**< pointer to time series */
   )
{
   assert(scip != NULL);
   assert(timeseries != NULL);

   BMSfreeMemoryArray(&(*timeseries)->name);

   SCIPfreeMemoryArray(scip, &(*timeseries)->vals);
   SCIPfreeMemoryArray(scip, &(*timeseries)->estimation);

   SCIPfreeMemory(scip, timeseries);

   *timeseries = NULL;
}

/** get current value of time series */
static
SCIP_Real timeSeriesGetValue(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   assert(timeseries != NULL);

   return timeseries->currentvalue;
}

/** get target value (which this time series reaches at the end of the solution process) */
static
SCIP_Real timeSeriesGetTargetValue(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   return timeseries->targetvalue;
}

/** get resolution of time series */
static
int timeSeriesGetResolution(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   return timeseries->resolution;
}

/** estimate tree size at which time series reaches target value */
static
SCIP_Real timeSeriesEstimate(
   TIMESERIES*           timeseries,         /**< time series */
   TREEDATA*             treedata            /**< tree data for fallback estimation */
   )
{
   SCIP_Real val;
   SCIP_Real targetval;
   SCIP_Real trend;
   SCIP_Real estimated;
   const SCIP_Real tolerance = 1e-6;

   /* if no observations have been made yet, return infinity */
   if( timeseries->nobs == 0L )
      return -1.0;

   val = timeSeriesGetValue(timeseries);
   targetval = timeSeriesGetTargetValue(timeseries);

   /* if the value has reached the target value already, return the number of observations */
   if( EPSZ(val - targetval, tolerance) )
      return treedata->nnodes;

   trend = doubleExpSmoothGetTrend(&timeseries->des);

   /* Get current value and trend. The linear trend estimation may point into the wrong direction
    * In this case, we use the fallback mechanism that we will need twice as many nodes.
    */
   if( (targetval > val && trend < tolerance) || (targetval < val && trend > -tolerance) )
   {
      return 2.0 * treedata->nvisited;
   }

   /* compute after how many additional steps the current trend reaches the target value; multiply by resolution */
   estimated = timeSeriesGetResolution(timeseries) * (timeseries->nvals + (targetval - val) / (SCIP_Real)trend);
   return timeseries->useleafts ? 2.0 * estimated - 1.0 : estimated;
}

/** update time series smoothened estimation */
static
void timeSeriesUpdateSmoothEstimation(
   TIMESERIES*           timeseries,         /**< time series */
   SCIP_Real             estimation          /**< estimation value */
   )
{
   if( timeseries->smoothestimation == SCIP_INVALID )/*lint !e777*/
      timeseries->smoothestimation = estimation;
   else
   {
      timeseries->smoothestimation *= (1.0 - SESCOEFF);
      timeseries->smoothestimation += SESCOEFF * estimation;
   }
}

/** get smooth estimation of time series */
static
SCIP_Real timeSeriesGetSmoothEstimation(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   return timeseries->smoothestimation;
}

/** resample to lower resolution */
static
void timeSeriesResample(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   DOUBLEEXPSMOOTH* des;
   int i;

   assert(timeseries->nvals % 2 == 0);

   des = &timeseries->des;
   doubleExpSmoothReset(des, timeseries->initialvalue);

   /* compress vals array to store only every second entry */
   for( i = 0; i < timeseries->nvals / 2; ++i )
   {
      timeseries->vals[i] = timeseries->vals[2 * i];
      timeseries->estimation[i] = timeseries->estimation[2 * i];
      doubleExpSmoothUpdate(des, timeseries->vals[i]);
      timeSeriesUpdateSmoothEstimation(timeseries, timeseries->estimation[i]);
   }

   timeseries->resolution *= 2;
   timeseries->nvals = timeseries->nvals / 2;
}

/** update time series */
static
SCIP_RETCODE timeSeriesUpdate(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMESERIES*           timeseries,         /**< time series */
   TREEDATA*             treedata,           /**< tree data */
   SCIP_Bool             isleaf              /**< are we at a leaf node? */
   )
{
   SCIP_Real value;

   assert(scip != NULL);
   assert(timeseries != NULL);
   assert(treedata != NULL);

   /* call update callback */
   assert(timeseries->timeseriesupdate != NULL);
   SCIP_CALL( timeseries->timeseriesupdate(scip, timeseries, treedata, &value) );

   /* store the value as current value */
   timeseries->currentvalue = value;

   if( timeseries->useleafts && ! isleaf )
      return SCIP_OKAY;

   timeseries->nobs++;

   /* if this is a leaf that matches the time series resolution, store the value */
   if( timeseries->nobs % timeseries->resolution == 0 )
   {
      int tspos;
      SCIP_Real estimate;

      assert(timeseries->nvals < timeseries->valssize);
      tspos = timeseries->nvals++;
      timeseries->vals[tspos] = value;
      doubleExpSmoothUpdate(&timeseries->des, value);
      estimate = timeSeriesEstimate(timeseries, treedata);
      timeseries->estimation[tspos] = estimate;
      timeSeriesUpdateSmoothEstimation(timeseries, estimate);
   }

   /* if the time series has reached its capacity, resample and increase the resolution */
   if( timeseries->nvals == timeseries->valssize )
      timeSeriesResample(timeseries);

   return SCIP_OKAY;
}

/** get name of time series */
static
char* timeSeriesGetName(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   return timeseries->name;
}

/** reset all time series */
static
void resetTimeSeries(
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   TIMESERIES** tss = eventhdlrdata->timeseries;
   int t;

   /* loop over time series and reset them */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      assert(tss[t] != NULL);
      timeSeriesReset(tss[t]);

      tss[t]->useleafts = eventhdlrdata->useleafts;
   }
}

/*
 * Callback methods of event handler
 */

/** free all time series */
static
void freeTimeSeries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   TIMESERIES** tss = eventhdlrdata->timeseries;
   int t;

   /* loop over time series and reset them */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      assert(tss[t] != NULL);
      timeSeriesFree(scip, &tss[t]);
   }
}

/** get ensemble tree size estimation as a combination of the individual time series estimations
 *
 *  the coefficients have been computed based on a nonlinear fit on a broad set of publicly available
 *  MIP instances; please refer to the publication at the top of this file for further details.
 */
static
SCIP_Real getEnsembleEstimation(
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   TREEDATA* treedata;
   SCIP_Real* coeffs;
   SCIP_Real estim;
   int t;

   TSPOS tsposs[] = {
      TSPOS_GAP,
      TSPOS_TREEWEIGHT,
      TSPOS_LFREQ,
      TSPOS_SSG,
      TSPOS_OPEN
   };

   /* coefficients for the early stage (tree weight <= 0.3) */
   SCIP_Real coeffs_early[] = {
      0.002, /* gap */
      0.381, /* tree weight */
      0.469, /* leaf-frequency */
      0.292, /* SSG */
      0.004  /* open-nodes */
   };

   /* coefficients for the intermediate stage (0.3 < tree weight <= 0.6) */
   SCIP_Real coeffs_intermediate[] = {
      0.011, /* gap */
      0.193, /* tree weight */
      0.351, /* leaf-frequency */
      0.012, /* SSG */
      0.051  /* open-nodes */
   };

   /* coefficients for the late stage (tree weight > 0.6) */
   SCIP_Real coeffs_late[] = {
      0.000, /* gap */
      0.033, /* tree weight */
      0.282, /* leaf-frequency */
      0.003, /* SSG */
      0.024  /* open-nodes */
   };

   assert(eventhdlrdata != NULL);
   treedata = eventhdlrdata->treedata;

   /* assign coeffs based on stage */
   if( treedata->weight <= 0.3 )
   {
      estim = 0.0;
      coeffs = coeffs_early;
      /* ensure that coeffs and time series are still aligned */
      assert(sizeof(coeffs_early)/sizeof(SCIP_Real) == NTIMESERIES); /*lint !e506*/
   }
   else if( treedata->weight <= 0.6 )
   {
      coeffs = coeffs_intermediate;
      /* ensure that coeffs and time series are still aligned */
      assert(sizeof(coeffs_intermediate)/sizeof(SCIP_Real) == NTIMESERIES); /*lint !e506*/

      /* initialize by intermediate WBE coefficient */
      estim = 0.156 * treeDataGetWbe(treedata);
   }
   else
   {
      coeffs = coeffs_late;
      /* ensure that coeffs and time series are still aligned */
      assert(sizeof(coeffs_late)/sizeof(SCIP_Real) == NTIMESERIES); /*lint !e506*/

      /* initialize by late WBE coefficient */
      estim = 0.579 * treeDataGetWbe(treedata);
   }

   /* combine estimation using the stage-dependent coefficients */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      SCIP_Real testim;
      TSPOS tspos = tsposs[t];
      testim = timeSeriesEstimate(eventhdlrdata->timeseries[tspos], treedata);

      if( testim < 0.0 )
         testim = treedata->nnodes;

      estim += coeffs[t] * testim;
   }

   if( estim < treedata->nnodes )
      return (SCIP_Real)treedata->nnodes;
   else
      return estim;
}

/** get approximation of search tree completion depending on the selected method */
static
SCIP_RETCODE getSearchCompletion(
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_Real*            completed           /**< pointer to store the search tree completion */
   )
{
   SCIP_Real values[9];
   TREEDATA* treedata;
   char completiontype;

   assert(eventhdlrdata != NULL);
   treedata = eventhdlrdata->treedata;
   completiontype = eventhdlrdata->completiontypeparam;

   /* infer automatic completion type
    *
    *  use regression forest if available,
    *  or
    *  use monotone regression if both SSG and tree weight are meaningful;
    *  or
    *  use tree weight or SSG, depending which one is available,
    *  or
    *  use gap, which is always available
    */
   if( completiontype == COMPLETIONTYPE_AUTO )
   {
      SCIP_Bool useweight = eventhdlrdata->treeisbinary;
      SCIP_Bool usessg = treedata->ssg->pblastsplit != SSG_STARTPRIMBOUND;/*lint !e777*/

      if( eventhdlrdata->regforest != NULL )
         completiontype = COMPLETIONTYPE_REGFOREST;
      else if( useweight && usessg )
         completiontype = COMPLETIONTYPE_MONOREG;
      else if( useweight )
         completiontype = COMPLETIONTYPE_TREEWEIGHT;
      else if( usessg )
         completiontype = COMPLETIONTYPE_SSG;
      else
         completiontype = COMPLETIONTYPE_GAP;
   }

   /* compute the search tree completion based on the selected method */
   switch (completiontype)
   {
   /* use regression forest */
   case COMPLETIONTYPE_REGFOREST:
      values[0] = timeSeriesGetValue(eventhdlrdata->timeseries[TSPOS_TREEWEIGHT]);
      values[1] = doubleExpSmoothGetTrend(&eventhdlrdata->timeseries[TSPOS_TREEWEIGHT]->des);
      values[2] = timeSeriesGetValue(eventhdlrdata->timeseries[TSPOS_SSG]);
      values[3] = doubleExpSmoothGetTrend(&eventhdlrdata->timeseries[TSPOS_SSG]->des);
      values[4] = timeSeriesGetValue(eventhdlrdata->timeseries[TSPOS_LFREQ]);
      values[5] = doubleExpSmoothGetTrend(&eventhdlrdata->timeseries[TSPOS_LFREQ]->des);
      values[6] = timeSeriesGetValue(eventhdlrdata->timeseries[TSPOS_GAP]);
      values[7] = doubleExpSmoothGetTrend(&eventhdlrdata->timeseries[TSPOS_GAP]->des);
      values[8] = doubleExpSmoothGetTrend(&eventhdlrdata->timeseries[TSPOS_OPEN]->des) < 0 ? 1.0 : 0.0;

      *completed = SCIPregForestPredict(eventhdlrdata->regforest, values);
      break;

   /* interpolate between ssg and tree weight */
   case COMPLETIONTYPE_MONOREG:
      *completed = eventhdlrdata->coefmonoweight * (SCIP_Real)treedata->weight +
         eventhdlrdata->coefmonossg * (1.0 - treedata->ssg->value);
      break;

   case COMPLETIONTYPE_TREEWEIGHT:
      *completed = (SCIP_Real)treedata->weight;
      break;

   case COMPLETIONTYPE_GAP:
      *completed = timeSeriesGetValue(eventhdlrdata->timeseries[TSPOS_GAP]); /* gap is stored as 1 - gap */
      break;

   case COMPLETIONTYPE_SSG:
      *completed = 1.0 - treedata->ssg->value; /* ssg is decreasing */
      break;

   default:
      SCIPerrorMessage("Unsupported completion type '%c'\n", completiontype);
      SCIPABORT();
      return SCIP_PARAMETERWRONGVAL;
   }
   return SCIP_OKAY;
}

/** tree size estimation based on search tree completion */
static
SCIP_RETCODE getEstimCompletion(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_Real*            estim               /**< pointer to store the estimation value */
   )
{
   SCIP_Real completed;

   *estim = -1.0;

   SCIP_CALL( getSearchCompletion(eventhdlrdata, &completed) );

   completed = MIN(completed, 1.0);

   if( completed > 0.0 )
      *estim = SCIPgetNNodes(scip) / completed;

   return SCIP_OKAY;
}

/** update callback at nodes */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateGap)
{ /*lint --e{715}*/
   SCIP_Real primalbound;
   SCIP_Real dualbound;

   assert(scip != NULL);
   assert(ts != NULL);
   assert(value != NULL);

   /* avoid to call SCIPgetDualbound during a restart where the queue is simply emptied */
   if( SCIPisInRestart(scip) )
   {
      *value = timeSeriesGetValue(ts);

      return SCIP_OKAY;
   }

   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);
   if( SCIPisInfinity(scip, REALABS(primalbound)) || SCIPisInfinity(scip, REALABS(dualbound)) )
      *value = 0;
   else if( SCIPisEQ(scip, primalbound, dualbound) )
      *value = 1.0;
   else
   {
      SCIP_Real abspb;
      SCIP_Real absdb;

      abspb = REALABS(primalbound);
      absdb = REALABS(dualbound);
      *value = 1.0 - REALABS(primalbound - dualbound)/MAX(abspb, absdb);
   }

   /* using this max, we set the closed gap to 0 in the case where the primal and dual bound differ in their sign */
   *value = MAX(*value, 0.0);

   return SCIP_OKAY;
}

/** update callback at nodes */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateTreeWeight)
{ /*lint --e{715}*/
   *value = (SCIP_Real)treedata->weight;

   return SCIP_OKAY;
}

/** update callback at nodes */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateLeafFreq)
{ /*lint --e{715}*/
   if( treedata->nvisited == 0 )
      *value = -0.5;
   else
      *value = (treedata->nleaves - 0.5)/(SCIP_Real)treedata->nvisited;

   return SCIP_OKAY;
}

/** update callback at nodes */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateSsg)
{ /*lint --e{715}*/
   if( treedata->nvisited == 0 )
      *value = 1.0;
   else
      *value = treedata->ssg->value;

   return SCIP_OKAY;
}

/** update callback at nodes */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateOpenNodes)
{ /*lint --e{715}*/
   if( treedata->nvisited == 0 )
      *value = 0.0;
   else
      *value = (SCIP_Real)treedata->nopen;

   return SCIP_OKAY;
}

/** include time series to forecast into event handler */
static
SCIP_RETCODE includeTimeseries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   /* include gap time series */
   SCIP_CALL( timeSeriesCreate(scip, &eventhdlrdata->timeseries[TSPOS_GAP], "gap", 1.0, 0.0,
            DES_ALPHA_GAP, DES_BETA_GAP, timeseriesUpdateGap) );

   /* include tree weight time series */
   SCIP_CALL( timeSeriesCreate(scip, &eventhdlrdata->timeseries[TSPOS_TREEWEIGHT], "tree-weight", 1.0, 0.0,
            DES_ALPHA_TREEWEIGHT, DES_BETA_TREEWEIGHT, timeseriesUpdateTreeWeight) );

   /* include leaf time series */
   SCIP_CALL( timeSeriesCreate(scip, &eventhdlrdata->timeseries[TSPOS_LFREQ], "leaf-frequency", 0.5, -0.5,
            DES_ALPHA_LEAFFREQUENCY, DES_BETA_LEAFFREQUENCY, timeseriesUpdateLeafFreq) );

   /* include SSG time series */
   SCIP_CALL( timeSeriesCreate(scip, &eventhdlrdata->timeseries[TSPOS_SSG], "ssg", 0.0, 1.0,
            DES_ALPHA_SSG, DES_BETA_SSG, timeseriesUpdateSsg) );

   /* include open nodes time series */
   SCIP_CALL( timeSeriesCreate(scip, &eventhdlrdata->timeseries[TSPOS_OPEN], "open-nodes", 0.0, 0.0,
            DES_ALPHA_OPENNODES, DES_BETA_OPENNODES, timeseriesUpdateOpenNodes) );

   return SCIP_OKAY;
}

/** get restartpolicy based on the value of the restart parameter */
static
RESTARTPOLICY getRestartPolicy(
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   switch (eventhdlrdata->restartpolicyparam)
   {
   case RESTARTPOLICY_CHAR_ALWAYS:
      return RESTARTPOLICY_ALWAYS;
   case RESTARTPOLICY_CHAR_NEVER:
      return RESTARTPOLICY_NEVER;
   case RESTARTPOLICY_CHAR_COMPLETION:
      return RESTARTPOLICY_COMPLETION;
   case RESTARTPOLICY_CHAR_ESTIMATION:
      return RESTARTPOLICY_ESTIMATION;
   default:
      SCIPerrorMessage("Unknown restart policy %c\n", eventhdlrdata->restartpolicyparam);
      break;
   }

   return RESTARTPOLICY_NEVER;
}

/** check if a restart is applicable considering limit and threshold user parameters */
static
SCIP_Bool isRestartApplicable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Longint nnodes;

   /* check whether to apply restarts when there are active pricers available */
   if( SCIPgetNActivePricers(scip) > 0 && ! eventhdlrdata->restartactpricers )
      return FALSE;

   /* check whether to apply a restart when nonlinear constraints are present */
   if( SCIPisNLPConstructed(scip) && ! eventhdlrdata->restartnonlinear )
      return FALSE;

   /* check if max number of restarts has been reached */
   if( eventhdlrdata->restartlimit != -1 && eventhdlrdata->nrestartsperformed >= eventhdlrdata->restartlimit )
      return FALSE;

   /* check if number of nodes exceeds the minimum number of nodes */
   if( eventhdlrdata->countonlyleaves )
      nnodes = eventhdlrdata->treedata->nleaves;
   else
      nnodes = eventhdlrdata->treedata->nvisited;

   if( nnodes < eventhdlrdata->minnodes )
      return FALSE;

   return TRUE;
}

/** should a restart be applied based on the value of the selected completion method? */   /*lint --e{715}*/
static
SCIP_Bool shouldApplyRestartCompletion(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{  /*lint --e{715}*/
   SCIP_Real completion;

   SCIP_CALL_ABORT( getSearchCompletion(eventhdlrdata, &completion) );

   /* if the estimation exceeds the current number of nodes by a dramatic factor, restart */
   if( completion < 1.0 / eventhdlrdata->restartfactor )
   {
      SCIPstatistic(
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "Completion %.5f less than restart threshold %.5f\n",
            completion, 1.0 / eventhdlrdata->restartfactor);
      )

      return TRUE;
   }

   return FALSE;
}

/** should a restart be applied based on the value of the selected completion method? */
static
SCIP_Bool shouldApplyRestartEstimation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Real estimation;

   estimation = SCIPgetTreesizeEstimation(scip);

   if( estimation < 0.0 )
   {
      SCIPstatistic(
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "Estimation %g is still unavailable\n",
            estimation);
      )

      return TRUE;
   }

   /* if the estimation exceeds the current number of nodes by a dramatic factor, restart */
   if( estimation > eventhdlrdata->treedata->nnodes * eventhdlrdata->restartfactor )
   {
      SCIPstatistic(
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "Estimation %g exceeds number of estimation tree nodes %" SCIP_LONGINT_FORMAT " by a factor of %.1f\n",
            estimation, eventhdlrdata->treedata->nnodes, estimation / eventhdlrdata->treedata->nnodes);
      )

      return TRUE;
   }

   return FALSE;
}

/** check if a restart should be performed based on the given restart policy */
static
SCIP_Bool shouldApplyRestart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Bool applyrestart = FALSE;

   switch (getRestartPolicy(eventhdlrdata))
   {
   case RESTARTPOLICY_ALWAYS:
      applyrestart = TRUE;
      break;
   case RESTARTPOLICY_NEVER:
      applyrestart = FALSE;
      break;
   case RESTARTPOLICY_COMPLETION:
      applyrestart = shouldApplyRestartCompletion(scip, eventhdlrdata);
      break;
   case RESTARTPOLICY_ESTIMATION:
      applyrestart = shouldApplyRestartEstimation(scip, eventhdlrdata);
      break;
   default:
      break;
   }

   return applyrestart;
}

/** update all time series */
static
SCIP_RETCODE updateTimeseries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   TREEDATA*             treedata,           /**< tree data */
   SCIP_Bool             isleaf              /**< are we at a leaf node? */
   )
{
   TIMESERIES** tss = eventhdlrdata->timeseries;
   int t;

   /* loop over time series */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      assert(tss[t] != NULL);
      SCIP_CALL( timeSeriesUpdate(scip, tss[t], treedata, isleaf) );

#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(scip,
         "Update of time series '%s', current value %.4f (%" SCIP_LONGINT_FORMAT " observations)\n",
         timeSeriesGetName(tss[t]), timeSeriesGetValue(tss[t]), tss[t]->nobs);
#endif
   }

   return SCIP_OKAY;
}

/** print a treesize estimation report into the string buffer */
static
char* printReport(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   char*                 strbuf,             /**< string buffer */
   int                   reportnum           /**< report number, or 0 to omit number */
   )
{
   TREEDATA* treedata = eventhdlrdata->treedata;
   char* ptr = strbuf;
   SCIP_Real completed;
   SCIP_Real wbeestim;
   char wbeestimstr[SCIP_MAXSTRLEN];
   int t;

   /* print report number */
   if( reportnum > 0 )
      ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN, "Report %d\nTime Elapsed: %.2f\n", reportnum, SCIPgetSolvingTime(scip));

   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
            "Estim. Tree Size   :%11" SCIP_LONGINT_FORMAT "\n",
            (SCIP_Longint)SCIPgetTreesizeEstimation(scip));

   SCIP_CALL_ABORT( getSearchCompletion(eventhdlrdata, &completed) );

   completed = MIN(1.0, completed);
   completed = MAX(0.0, completed);

   /* print tree data */
   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN,
         "%-19s: %" SCIP_LONGINT_FORMAT " nodes ("
         "%" SCIP_LONGINT_FORMAT " visited, "
         "%" SCIP_LONGINT_FORMAT " internal, "
         "%" SCIP_LONGINT_FORMAT " leaves, "
         "%" SCIP_LONGINT_FORMAT " open), "
         "weight: %.4Lf completed %.4f\n",
         "Estimation Tree",
         treedata->nnodes,
         treedata->nvisited,
         treedata->ninner,
         treedata->nleaves,
         treedata->nopen,
         treedata->weight,
         completed
         );

   /* print estimations */
   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN, "Estimations        : %10s %10s %10s %10s %10s",
            "estim", "value", "trend", "resolution", "smooth");
   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN, "\n");

   wbeestim = treeDataGetWbe(eventhdlrdata->treedata);
   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN, "  wbe              : %10s %10s %10s %10s %10s\n",
            real2String(wbeestim, wbeestimstr, 0), "-", "-", "-", "-");

   ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN, "  tree-profile     : %10.0f %10s %10s %10s %10s\n",
            predictTotalSizeTreeProfile(scip, eventhdlrdata->treeprofile, eventhdlrdata->treeprofile_minnodesperdepth),
            "-", "-", "-", "-");

   /* print time series forecasts */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      SCIP_Real trend;
      SCIP_Real smoothestim;
      TIMESERIES* ts = eventhdlrdata->timeseries[t];
      char trendstr[SCIP_MAXSTRLEN];
      char smoothestimstr[SCIP_MAXSTRLEN];

      trend = doubleExpSmoothGetTrend(&ts->des);
      smoothestim = timeSeriesGetSmoothEstimation(ts);

      ptr += SCIPsnprintf(ptr, SCIP_MAXSTRLEN, "  %-17s: %10.0f %10.5f %10s %10d %10s\n",
            timeSeriesGetName(ts),
            timeSeriesEstimate(ts, eventhdlrdata->treedata),
            timeSeriesGetValue(ts),
            real2String(trend, trendstr, 5),
            timeSeriesGetResolution(ts),
            real2String(smoothestim, smoothestimstr, 0));
   }

   if( reportnum > 0 )
      (void) SCIPsnprintf(ptr, SCIP_MAXSTRLEN, "End of Report %d\n", reportnum);

   return strbuf;
}


/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyEstim)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeEventHdlrEstim(scip) );

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeEstim)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   freeTreeData(scip, &eventhdlrdata->treedata);

   freeTimeSeries(scip, eventhdlrdata);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitEstim)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* test if user specified a regression forest */
   if( 0 != strncmp(eventhdlrdata->regforestfilename, DEFAULT_REGFORESTFILENAME, strlen(DEFAULT_REGFORESTFILENAME)) )
   {
      SCIP_CALL( SCIPregForestFromFile(&eventhdlrdata->regforest, eventhdlrdata->regforestfilename) );
   }

   eventhdlrdata->lastrestartrun = 0;
   eventhdlrdata->nrestartsperformed = 0;

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitEstim)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPregForestFree(&eventhdlrdata->regforest);

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolEstim)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->restarthitcounter = 0;
   eventhdlrdata->weightlastreport = 0.0;
   eventhdlrdata->nreports = 0;

   /* reset tree data */
   SCIP_CALL( resetTreeData(scip, eventhdlrdata->treedata) );

   resetTimeSeries(eventhdlrdata);

   SCIP_CALL( SCIPcatchEvent(scip, EVENTTYPE_ESTIM, eventhdlr, NULL, NULL) );

   if( eventhdlrdata->treeprofile_enabled )
   {
      SCIP_CALL( createTreeProfile(scip, &eventhdlrdata->treeprofile) );
   }

   eventhdlrdata->treeisbinary = TRUE;

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolEstim)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->treeprofile != NULL )
      freeTreeProfile(scip, &eventhdlrdata->treeprofile);

   SCIP_CALL( SCIPdropEvent(scip, EVENTTYPE_ESTIM, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecEstim)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_Bool isleaf;
   SCIP_EVENTTYPE eventtype;
   TREEDATA* treedata;
   char strbuf[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   eventtype = SCIPeventGetType(event);
   treedata = eventhdlrdata->treedata;

   /* actual leaf nodes for our tree data are children/siblings/leaves or the focus node itself (deadend)
    * if it has not been branched on
    */
   isleaf = (eventtype == SCIP_EVENTTYPE_NODEDELETE) &&
      (SCIPnodeGetType(SCIPeventGetNode(event)) == SCIP_NODETYPE_CHILD ||
         SCIPnodeGetType(SCIPeventGetNode(event)) == SCIP_NODETYPE_SIBLING ||
         SCIPnodeGetType(SCIPeventGetNode(event)) == SCIP_NODETYPE_LEAF ||
         (SCIPnodeGetType(SCIPeventGetNode(event)) == SCIP_NODETYPE_DEADEND && !SCIPwasNodeLastBranchParent(scip, SCIPeventGetNode(event))));

   if( eventtype == SCIP_EVENTTYPE_NODEBRANCHED || isleaf )
   {
      SCIP_NODE* eventnode;
      int nchildren = 0;

      if( eventtype == SCIP_EVENTTYPE_NODEBRANCHED )
      {
         nchildren = SCIPgetNChildren(scip);

         /* update whether the tree is still binary */
         if( nchildren != 2 )
            eventhdlrdata->treeisbinary = FALSE;
      }

      eventnode = SCIPeventGetNode(event);
      SCIP_CALL( updateTreeData(scip, treedata, eventnode, nchildren) );
      SCIP_CALL( updateTreeProfile(scip, eventhdlrdata->treeprofile, eventnode) );

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "%s\n", treeDataPrint(treedata, strbuf));
#endif

      SCIP_CALL( updateTimeseries(scip, eventhdlrdata, treedata, nchildren == 0) );

      /* should a new report be printed? */
      if( eventhdlrdata->reportfreq >= 0 && SCIPgetStatus(scip) == SCIP_STATUS_UNKNOWN &&
         (eventhdlrdata->reportfreq == 0
         || treedata->weight >= eventhdlrdata->weightlastreport + 1.0 / (SCIP_Real)eventhdlrdata->reportfreq) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s\n", printReport(scip, eventhdlrdata, strbuf, ++eventhdlrdata->nreports));

         if( eventhdlrdata->reportfreq > 0 )
            eventhdlrdata->weightlastreport = 1 / (SCIP_Real)eventhdlrdata->reportfreq * SCIPfloor(scip, ((SCIP_Real)treedata->weight * eventhdlrdata->reportfreq));
         else
            eventhdlrdata->weightlastreport = (SCIP_Real)treedata->weight;
      }
   }

   /* if nodes have been pruned, things are progressing, don't restart right now */
   if( isleaf )
      return SCIP_OKAY;

   /* check if all conditions are met such that the event handler should run */
   if( ! isRestartApplicable(scip, eventhdlrdata) )
      return SCIP_OKAY;

   /* test if a restart should be applied */
   if( shouldApplyRestart(scip, eventhdlrdata) )
   {
      eventhdlrdata->restarthitcounter++;

      if( eventhdlrdata->restarthitcounter >= eventhdlrdata->hitcounterlim )
      {
         /* safe that we triggered a restart at this run */
         if( SCIPgetNRuns(scip) > eventhdlrdata->lastrestartrun )
         {
            eventhdlrdata->nrestartsperformed++;

            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "Restart triggered after %d consecutive estimations that the remaining tree will be large\n",
               eventhdlrdata->restarthitcounter);
         }

         eventhdlrdata->lastrestartrun = SCIPgetNRuns(scip);

         SCIP_CALL( SCIPrestartSolve(scip) );
      }
   }
   else
   {
      eventhdlrdata->restarthitcounter = 0;
   }

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputEstim)
{  /*lint --e{715}*/
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   char strbuf[SCIP_MAXSTRLEN];

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->showstats )
      SCIPinfoMessage(scip, file, "%s", printReport(scip, eventhdlrdata, strbuf, 0));

   return SCIP_OKAY;
}

/** output method of search tree completion display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputCompleted)
{  /*lint --e{715}*/
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   TREEDATA* treedata;
   SCIP_Real completed;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME) == 0);
   assert(scip != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   treedata = eventhdlrdata->treedata;

   SCIP_CALL( getSearchCompletion(eventhdlrdata, &completed) );

   completed = MIN(completed, 1.0);

   if( treedata->weight >= 0.005 && completed > 0 )
      SCIPinfoMessage(scip, file, "%7.2f%%", 100.0 * completed);
   else
      SCIPinfoMessage(scip, file, " unknown");

   return SCIP_OKAY;
}

/** creates event handler for tree size estimation */
SCIP_RETCODE SCIPincludeEventHdlrEstim(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RETCODE retcode;
   SCIP_EVENTHDLRDATA* eventhdlrdata = NULL;
   SCIP_EVENTHDLR* eventhdlr = NULL;

   /* create estim event handler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   BMSclearMemory(eventhdlrdata);

   SCIP_CALL_TERMINATE( retcode, createTreeData(scip, &eventhdlrdata->treedata), TERMINATE );

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecEstim, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyEstim) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeEstim) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitEstim) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitEstim) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolEstim) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolEstim) );

   /* add estimation event handler parameters */
   SCIP_CALL( SCIPaddCharParam(scip, "estimation/restarts/restartpolicy", "restart policy: (a)lways, (c)ompletion, (e)stimation, (n)ever",
         &eventhdlrdata->restartpolicyparam, FALSE, DEFAULT_RESTARTPOLICY, "acen", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "estimation/method",
            "tree size estimation method: (c)ompletion, (e)nsemble, "
            "time series forecasts on either (g)ap, (l)eaf frequency, (o)open nodes, tree (w)eight, (s)sg, "
            "or (t)ree profile or w(b)e",
         &eventhdlrdata->estimmethod, FALSE, DEFAULT_ESTIMMETHOD, ESTIMMETHODS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "estimation/restarts/restartlimit", "restart limit",
         &eventhdlrdata->restartlimit, FALSE, DEFAULT_RESTARTLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "estimation/restarts/minnodes", "minimum number of nodes before restart",
         &eventhdlrdata->minnodes, FALSE, DEFAULT_MINNODES, -1L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "estimation/restarts/countonlyleaves", "should only leaves count for the minnodes parameter?",
         &eventhdlrdata->countonlyleaves, DEFAULT_COUNTONLYLEAVES, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "estimation/restarts/restartfactor",
         "factor by which the estimated number of nodes should exceed the current number of nodes",
         &eventhdlrdata->restartfactor, FALSE, DEFAULT_RESTARTFACTOR, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "estimation/restarts/restartnonlinear",
         "whether to apply a restart when nonlinear constraints are present",
         &eventhdlrdata->restartnonlinear, FALSE, DEFAULT_RESTARTNONLINEAR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "estimation/restarts/restartactpricers",
         "whether to apply a restart when active pricers are used",
         &eventhdlrdata->restartactpricers, FALSE, DEFAULT_RESTARTACTPRICERS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "estimation/coefmonoweight",
         "coefficient of tree weight in monotone approximation of search completion",
         &eventhdlrdata->coefmonoweight, FALSE, DEFAULT_COEFMONOWEIGHT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "estimation/coefmonossg",
         "coefficient of 1 - SSG in monotone approximation of search completion",
         &eventhdlrdata->coefmonossg, FALSE, DEFAULT_COEFMONOSSG, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "estimation/restarts/hitcounterlim", "limit on the number of successive samples to really trigger a restart",
         &eventhdlrdata->hitcounterlim, FALSE, DEFAULT_HITCOUNTERLIM, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "estimation/reportfreq",
         "report frequency on estimation: -1: never, 0:always, k >= 1: k times evenly during search",
         &eventhdlrdata->reportfreq, TRUE, DEFAULT_REPORTFREQ, -1, INT_MAX / 2, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "estimation/regforestfilename", "user regression forest in RFCSV format",
         &eventhdlrdata->regforestfilename, FALSE, DEFAULT_REGFORESTFILENAME, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "estimation/completiontype",
         "approximation of search tree completion: (a)uto, (g)ap, tree (w)eight, (m)onotone regression, (r)egression forest, (s)sg",
         &eventhdlrdata->completiontypeparam, FALSE, DEFAULT_COMPLETIONTYPE, "agmrsw", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "estimation/treeprofile/enabled",
         "should the event handler collect data?",
         &eventhdlrdata->treeprofile_enabled, FALSE, DEFAULT_TREEPROFILE_ENABLED, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "estimation/treeprofile/minnodesperdepth",
         "minimum average number of nodes at each depth before producing estimations",
         &eventhdlrdata->treeprofile_minnodesperdepth, FALSE, DEFAULT_TREEPROFILE_MINNODESPERDEPTH, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "estimation/useleafts",
         "use leaf nodes as basic observations for time series, or all nodes?",
         &eventhdlrdata->useleafts, TRUE, DEFAULT_USELEAFTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "estimation/showstats",
         "should statistics be shown at the end?",
         &eventhdlrdata->showstats, TRUE, DEFAULT_SHOWSTATS, NULL, NULL) );

   /* SSG parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "estimation/ssg/nmaxsubtrees",
         "the maximum number of individual SSG subtrees; -1: no limit",
         &eventhdlrdata->treedata->ssg->nmaxsubtrees, FALSE, DEFAULT_SSG_NMAXSUBTREES, -1, INT_MAX / 2, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "estimation/ssg/nminnodeslastsplit",
         "minimum number of nodes to process between two consecutive SSG splits",
         &eventhdlrdata->treedata->ssg->nminnodeslastsplit, FALSE, DEFAULT_SSG_NMINNODESLASTSPLIT, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME, TABLE_DESC, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputEstim,
         NULL, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   /* include time series into event handler */
   SCIP_CALL( includeTimeseries(scip, eventhdlrdata) );

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, SCIP_DISPSTATUS_AUTO,
         NULL, NULL, NULL, NULL, NULL, NULL, dispOutputCompleted,
         NULL, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

/* cppcheck-suppress unusedLabel */
TERMINATE:
   if( retcode != SCIP_OKAY )
   {
      freeTreeData(scip, &eventhdlrdata->treedata);
      SCIPfreeMemory(scip, &eventhdlrdata);
   }

   return retcode;
}

/** return an estimation of the final tree size */
SCIP_Real SCIPgetTreesizeEstimation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   TSPOS tspos = TSPOS_NONE;
   SCIP_Real estim;

   assert(scip != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( eventhdlr == NULL )
   {
      SCIPwarningMessage(scip, "SCIPgetTreesizeEstimation() called, but event handler " EVENTHDLR_NAME " is missing.\n");
      return -1.0;
   }

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   switch (eventhdlrdata->estimmethod)
   {
   case ESTIMMETHOD_COMPL:
      SCIP_CALL_ABORT( getEstimCompletion(scip, eventhdlrdata, &estim) );
      return estim;

   case ESTIMMETHOD_ENSMBL:
      return getEnsembleEstimation(eventhdlrdata);

   /* for the requested time series methods, we specify the array position */
   case ESTIMMETHOD_GAP:
      tspos = TSPOS_GAP;
      break;

   case ESTIMMETHOD_LFREQ:
      tspos = TSPOS_LFREQ;
      break;

   case ESTIMMETHOD_OPEN:
      tspos = TSPOS_OPEN;
      break;

   case ESTIMMETHOD_TREEWEIGHT:
      tspos = TSPOS_TREEWEIGHT;
      break;

   case ESTIMMETHOD_SSG:
      tspos = TSPOS_SSG;
      break;

   /* tree profile estimation */
   case ESTIMMETHOD_TPROF:
      return predictTotalSizeTreeProfile(scip, eventhdlrdata->treeprofile, eventhdlrdata->treeprofile_minnodesperdepth);

   /* Weighted backtrack estimation */
   case ESTIMMETHOD_WBE:
      return treeDataGetWbe(eventhdlrdata->treedata);

   default:
      SCIPerrorMessage("Unknown estimation '%c' method specified, should be one of [%s]\n",
         eventhdlrdata->estimmethod, ESTIMMETHODS);
      SCIPABORT();
      break;
   }

   assert(tspos != TSPOS_NONE);
   return (tspos == TSPOS_NONE ? -1.0 : timeSeriesEstimate(eventhdlrdata->timeseries[tspos], eventhdlrdata->treedata));
}
