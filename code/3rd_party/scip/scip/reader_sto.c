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

/**@file   reader_sto.c
 * @ingroup DEFPLUGINS_READER
 * @brief  STO file reader - the stochastic information of an instance in SMPS format
 * @author Stephen J. Maher
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/benders_default.h"
#include "scip/cons_linear.h"
#include "scip/pub_cons.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/pub_var.h"
#include "scip/reader_cor.h"
#include "scip/reader_sto.h"
#include "scip/reader_tim.h"
#include "scip/scip_cons.h"
#include "scip/scip_debug.h"
#include "scip/scipdefplugins.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include "scip/scip_var.h"
#include <stdlib.h>
#include <string.h>

#define READER_NAME             "storeader"
#define READER_DESC             "file reader for stochastic information of stochastic programs in the SMPS file format"
#define READER_EXTENSION        "sto"

#define DEFAULT_USEBENDERS            FALSE  /**< should Benders' decomposition be used for the stochastic program? */

/*
 * sto reader internal methods
 */

#define STO_MAX_LINELEN  1024
#define STO_MAX_NAMELEN   256

#define STO_DEFAULT_ARRAYSIZE          100
#define STO_DEFAULT_ENTRIESSIZE         20
#define STO_DEFAULT_BLOCKARRAYSIZE       5
#define STO_DEFAULT_CHILDRENSIZE         5

#define BLANK         ' '

typedef struct StoScenario STOSCENARIO;

/** STO reading data */
struct SCIP_ReaderData
{
   SCIP_Bool             usebenders;
   STOSCENARIO*          scenariotree;       /**< the multi stage scenario tree */
   int                   numscenarios;       /**< the total number of scenarios in the scenario tree */
};


struct StoScenario
{
   SCIP*                 scip;               /**< the SCIP instance for the scenario. Used for benders. */
   SCIP**                subproblems;        /**< the SCIP instances for the subproblems */
   STOSCENARIO*          parent;             /**< parent scenario. */
   STOSCENARIO**         children;           /**< children scenarios. */
   int                   nchildren;          /**< the number of children scenarios. */
   int                   childrensize;       /**< the size of the children array. */
   int                   nsubproblems;       /**< the number of subproblems */
   int                   stagenum;           /**< the number of the stage */
   int                   scenarionum;        /**< the scenario number of this stage */
   const char*           stagename;          /**< the stage name */
   const char*           name;               /**< the scenario name. */
   SCIP_Real             probability;        /**< the probability for this scenario. */
   SCIP_Real             lowerbound;         /**< the lower bound for this scenario */
   /* the following describes the modifications to the constraint matrix and rhs for each scenario. */
   const char**          rownames;           /**< the names of the rows with a changed value. */
   const char**          colnames;           /**< the names of the columns with a changed value. */
   SCIP_Real*            values;             /**< the values for the given row/column pair. */
   int                   nentries;           /**< the number of row/column pairs */
   int                   entriessize;        /**< the size of the row/colum arrays */
};



/** enum containing all sto sections */
enum StoSection
{
   STO_STOCH,
   STO_SCENARIOS,
   STO_BLOCKS,
   STO_INDEP,
   STO_ENDATA
};
typedef enum StoSection STOSECTION;

/** enum containing the types of stochastic information */
enum StoStochInfo
{
   STO_STOCHINFO_NONE         = -1,
   STO_STOCHINFO_DISCRETE     = 0,
   STO_STOCHINFO_UNIFORM      = 1,
   STO_STOCHINFO_NORMAL       = 2,
   STO_STOCHINFO_SUB          = 3,
   STO_STOCHINFO_LINTR        = 4
};
typedef enum StoStochInfo STOSTOCHINFO;

/** sto input structure */
struct StoInput
{
   STOSECTION            section;
   STOSTOCHINFO          stochinfotype;
   SCIP_FILE*            fp;
   int                   lineno;
   SCIP_Bool             haserror;
   char                  buf[STO_MAX_LINELEN];
   const char*           f0;
   const char*           f1;
   const char*           f2;
   const char*           f3;
   const char*           f4;
   const char*           f5;
   const char*           f6;
   char                  probname[STO_MAX_NAMELEN];
   char                  stochtype[STO_MAX_NAMELEN];
};
typedef struct StoInput STOINPUT;

/** creates a scenario structure */
static
SCIP_RETCODE createScenarioData(
   SCIP*                 scip,               /**< SCIP data structure */
   STOSCENARIO**         scenariodata        /**< the scenario to be created */
   )
{
   assert(scip != NULL);

   SCIPdebugMessage("Creating scenario data.\n");

   SCIP_CALL( SCIPallocBlockMemory(scip, scenariodata) );

   (*scenariodata)->scip = NULL;
   (*scenariodata)->subproblems = NULL;
   (*scenariodata)->parent = NULL;
   (*scenariodata)->nchildren = 0;
   (*scenariodata)->childrensize = STO_DEFAULT_CHILDRENSIZE;
   (*scenariodata)->nsubproblems = 0;
   (*scenariodata)->stagenum = -1;
   (*scenariodata)->scenarionum = -1;
   (*scenariodata)->stagename = NULL;
   (*scenariodata)->name = NULL;
   (*scenariodata)->probability = 1.0;
   (*scenariodata)->lowerbound = -SCIPinfinity(scip);
   (*scenariodata)->nentries = 0;
   (*scenariodata)->entriessize = STO_DEFAULT_ENTRIESSIZE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*scenariodata)->children, (*scenariodata)->childrensize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*scenariodata)->rownames, (*scenariodata)->entriessize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*scenariodata)->colnames, (*scenariodata)->entriessize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*scenariodata)->values, (*scenariodata)->entriessize) );

   return SCIP_OKAY;
}

/** frees the memory used for the scenario tree */
static
SCIP_RETCODE freeScenarioTree(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO**         scenariotree        /**< the scenario tree */
   )
{
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(scenariotree != NULL);
   assert(*scenariotree != NULL);

   SCIPdebugMessage("Freeing scenario <%s> in stage <%s>\n", (*scenariotree)->name,
      (*scenariotree)->stagename);

   /* storing the number of children before starting the recursive freeing */
   nchildren = (*scenariotree)->nchildren;

   while( (*scenariotree)->nchildren > 0 )
   {
      SCIP_CALL( freeScenarioTree(scip, &(*scenariotree)->children[(*scenariotree)->nchildren - 1]) );
      (*scenariotree)->nchildren--;
   }

   for( i = (*scenariotree)->nentries - 1; i >= 0; i-- )
   {
      SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->colnames[i], strlen((*scenariotree)->colnames[i]) + 1);
      SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->rownames[i], strlen((*scenariotree)->rownames[i]) + 1);
   }

   SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->values, (*scenariotree)->entriessize);
   SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->colnames, (*scenariotree)->entriessize);
   SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->rownames, (*scenariotree)->entriessize);
   SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->children, (*scenariotree)->childrensize);

   SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->name, strlen((*scenariotree)->name) + 1);
   SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->stagename, strlen((*scenariotree)->stagename) + 1);

   /* freeing the subproblem SCIP instances */
   for( i = (*scenariotree)->nsubproblems - 1; i >= 0; i-- )
      SCIP_CALL( SCIPfree(&(*scenariotree)->subproblems[i]) );

   /* freeing the array that stores the subproblems */
   if( nchildren > 0 && (*scenariotree)->subproblems != NULL )
      SCIPfreeBlockMemoryArray(scip, &(*scenariotree)->subproblems, nchildren);

   SCIPfreeBlockMemory(scip, scenariotree);

   return SCIP_OKAY;
}

/** sets the SCIP pointer to the scenario */
static
void setScenarioScip(
   STOSCENARIO*          scenario,           /**< the scenario */
   SCIP*                 scip                /**< the SCIP data structure */
   )
{
   assert(scenario != NULL);
   assert(scip != NULL);

   scenario->scip = scip;
}

/** returns the SCIP pointer to the scenario */
static
SCIP* getScenarioScip(
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scenario != NULL);

   return scenario->scip;
}

/** creates the subproblem array. This array will be the same size as the number of children */
static
SCIP_RETCODE createScenarioSubproblemArray(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scenario->subproblems, scenario->nchildren) );

   return SCIP_OKAY;
}

/** adds a scenario to the subproblem array */
static
void addScenarioSubproblem(
   STOSCENARIO*          scenario,           /**< the scenario */
   SCIP*                 subproblem          /**< the subproblems data structure */
   )
{
   assert(scenario != NULL);
   assert(subproblem != NULL);

   assert(scenario->nsubproblems + 1 <= scenario->nchildren);

   scenario->subproblems[scenario->nsubproblems] = subproblem;
   scenario->nsubproblems++;
}

/** returns the subproblem array for the scenario */
static
SCIP** getScenarioSubproblemArray(
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scenario != NULL);

   return scenario->subproblems;
}

/** returns the number of children for a given scenario */
static
int getScenarioNChildren(
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scenario != NULL);

   return scenario->nchildren;
}

/** returns a given child for a given scenario */
static
STOSCENARIO* getScenarioChild(
   STOSCENARIO*          scenario,           /**< the scenario */
   int                   childnum            /**< the number of the desired child */
   )
{
   assert(scenario != NULL);
   assert(childnum >= 0 && childnum < scenario->nchildren);

   return scenario->children[childnum];
}

/** returns the parent of a scenario */
static
STOSCENARIO* getScenarioParent(
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scenario != NULL);

   return scenario->parent;
}

/** sets the stage name */
static
SCIP_RETCODE setScenarioStageName(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the scenario */
   const char*           stagename           /**< the stage name */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &scenario->stagename, stagename, strlen(stagename) + 1) );

   return SCIP_OKAY;
}

/** returns the stage name */
static
const char* getScenarioStageName(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   return scenario->stagename;
}

/** sets the stage num */
static
SCIP_RETCODE setScenarioStageNum(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the scenario */
   int                   stagenum            /**< the stage num */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   scenario->stagenum = stagenum;

   return SCIP_OKAY;
}

/** returns the stage num */
static
int getScenarioStageNum(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   return scenario->stagenum;
}

/** sets the scenario name */
static
SCIP_RETCODE setScenarioName(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the scenario */
   const char*           name                /**< the scenario name */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &scenario->name, name, strlen(name) + 1) );

   return SCIP_OKAY;
}

/** returns the scenario name */
static
const char* getScenarioName(
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scenario != NULL);

   return scenario->name;
}

/** sets the scenario num */
static
SCIP_RETCODE setScenarioNum(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the scenario */
   int                   scenarionum         /**< the scenario num */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   scenario->scenarionum = scenarionum;

   return SCIP_OKAY;
}

/** returns the scenario num */
static
int getScenarioNum(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   return scenario->scenarionum;
}

/** sets the scenario probability */
static
SCIP_RETCODE setScenarioProbability(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the scenario */
   SCIP_Real             probability         /**< the scenario probability */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   scenario->probability = probability;

   return SCIP_OKAY;
}

/** returns the scenario probability */
static
SCIP_Real getScenarioProbability(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   return scenario->probability;
}

/** sets the scenario lowerbound */
static
SCIP_RETCODE setScenarioLowerbound(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the scenario */
   SCIP_Real             lowerbound          /**< the scenario lowerbound */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   scenario->lowerbound = lowerbound;

   return SCIP_OKAY;
}

/** returns the scenario lowerbound */
static
SCIP_Real getScenarioLowerbound(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   return scenario->lowerbound;
}

/** add scenario entry */
static
SCIP_RETCODE addScenarioEntry(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the scenario */
   const char*           rowname,            /**< the row name for the entry */
   const char*           colname,            /**< the col name for the entry */
   SCIP_Real             value               /**< the value for the entry */
   )
{
   assert(scip != NULL);
   assert(scenario != NULL);

   if( scenario->nentries + 1 > scenario->entriessize )
   {
      int newsize;
      newsize = SCIPcalcMemGrowSize(scip, scenario->nentries + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scenario->rownames, scenario->entriessize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scenario->colnames, scenario->entriessize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scenario->values, scenario->entriessize, newsize) );
      scenario->entriessize = newsize;
   }

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &scenario->rownames[scenario->nentries], rowname, strlen(rowname) + 1) );   /*lint !e866*/
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &scenario->colnames[scenario->nentries], colname, strlen(colname) + 1) );   /*lint !e866*/

   scenario->values[scenario->nentries] = value;
   scenario->nentries++;

   return SCIP_OKAY;
}

/** returns the number of entries for a scenario */
static
int getScenarioNEntries(
   STOSCENARIO*          scenario            /**< the scenario */
   )
{
   assert(scenario != NULL);

   return scenario->nentries;
}

/** returns an entry row for a scenario */
static
const char* getScenarioEntryRow(
   STOSCENARIO*          scenario,           /**< the scenario */
   int                   entry               /**< the entry number */
   )
{
   assert(scenario != NULL);
   assert(entry >= 0 && entry < scenario->nentries);

   return scenario->rownames[entry];
}

/** returns an entry column for a scenario */
static
const char* getScenarioEntryCol(
   STOSCENARIO*          scenario,           /**< the scenario */
   int                   entry               /**< the entry number */
   )
{
   assert(scenario != NULL);
   assert(entry >= 0 && entry < scenario->nentries);

   return scenario->colnames[entry];
}

/** returns an entry value for a scenario */
static
SCIP_Real getScenarioEntryValue(
   STOSCENARIO*          scenario,           /**< the scenario */
   int                   entry               /**< the entry number */
   )
{
   assert(scenario != NULL);
   assert(entry >= 0 && entry < scenario->nentries);

   return scenario->values[entry];
}

/** copies a scenario.
 * In the case of blocks, the scenarios must be combined
 */
static
SCIP_RETCODE copyScenario(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          sourcescenario,     /**< the source scenario */
   STOSCENARIO**         targetscenario,     /**< the target scenario */
   SCIP_Bool             copyname            /**< should the name be copied? */
   )
{
   SCIP_Real probability;
   SCIP_Real lowerbound;
   int i;

   assert(scip != NULL);
   assert(sourcescenario != NULL);
   assert(targetscenario != NULL);

   /* setting the stage name */
   if( copyname )
   {
      SCIP_CALL( setScenarioName(scip, (*targetscenario), sourcescenario->name) );
      SCIP_CALL( setScenarioStageName(scip, (*targetscenario), sourcescenario->stagename) );
      SCIP_CALL( setScenarioNum(scip, (*targetscenario), sourcescenario->scenarionum) );
      SCIP_CALL( setScenarioStageNum(scip, (*targetscenario), sourcescenario->stagenum) );
   }

   /* adding the entries from scenario 1 and 2 to the merged scenario */
   for( i = 0; i < sourcescenario->nentries; i++ )
      SCIP_CALL( addScenarioEntry(scip, (*targetscenario), sourcescenario->rownames[i], sourcescenario->colnames[i],
            sourcescenario->values[i]) );

   /* setting the scenario probability */
   probability = getScenarioProbability(scip, sourcescenario);
   SCIP_CALL( setScenarioProbability(scip, (*targetscenario), probability) );

   lowerbound = getScenarioLowerbound(scip, sourcescenario);
   SCIP_CALL( setScenarioLowerbound(scip, (*targetscenario), lowerbound) );

   return SCIP_OKAY;
}

/** merge scenarios.
 *  In the case of blocks, the scenarios must be combined
 */
static
SCIP_RETCODE mergeScenarios(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario1,          /**< the first scenario */
   STOSCENARIO**         mergedscenario      /**< the merged scenario */
   )
{
   SCIP_Real probability;
   int i;

   assert(scip != NULL);
   assert(scenario1 != NULL);
   assert(mergedscenario != NULL);

   /* adding the entries from scenario 1 and 2 to the merged scenario */
   for( i = 0; i < scenario1->nentries; i++ )
      SCIP_CALL( addScenarioEntry(scip, (*mergedscenario), scenario1->rownames[i], scenario1->colnames[i],
            scenario1->values[i]) );

   /* setting the scenario probability */
   probability = getScenarioProbability(scip, scenario1)*getScenarioProbability(scip, (*mergedscenario));
   SCIP_CALL( setScenarioProbability(scip, (*mergedscenario), probability) );

   return SCIP_OKAY;
}

/** adds a child to a given scenario */
static
SCIP_RETCODE scenarioAddChild(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO**         parent,             /**< the parent scenario */
   STOSCENARIO*          child               /**< the child scenario */
   )
{
   STOSCENARIO* scenario;

   assert(parent != NULL);
   assert((*parent) != NULL);
   assert(child != NULL);

   if( (*parent)->nchildren + 1 > (*parent)->childrensize )
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(*parent)->children, &(*parent)->childrensize,
            (*parent)->nchildren + 1) );

   SCIP_CALL( createScenarioData(scip, &scenario) );
   SCIP_CALL( copyScenario(scip, child, &scenario, TRUE) );
   scenario->parent = (*parent);

   (*parent)->children[(*parent)->nchildren] = scenario;
   (*parent)->nchildren++;

   return SCIP_OKAY;
}

/** recursively adds the scenarios to the reader data */
static
SCIP_RETCODE buildScenarioTree(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO**         scenariotree,       /**< the scenario tree */
   STOSCENARIO***        scenarios,          /**< the array of scenarios */
   int*                  numscenarios,       /**< the number of scenarios per stage */
   int                   numstages,          /**< the number of stages */
   int                   stage               /**< the number of the stage. Also the depth of the tree */
   )
{
   int stageindex;
   int i;

   assert(scip != NULL);
   assert(scenariotree != NULL);
   assert(stage >= 0 && stage < numstages);

   /* finding the scenarios for this stage */
   for( i = 0; i < numstages; i++ )
   {
      if( strcmp(getScenarioStageName(scip, scenarios[i][0]), SCIPtimGetStageName(scip, stage + 1)) == 0 )
         break;
   }
   assert(i < numstages);

   stageindex = i;

   /* adds each scenario to the scenario tree */
   for( i = 0; i < numscenarios[stageindex]; i++ )
   {
      /* adding child to the scenario tree */
      SCIP_CALL( scenarioAddChild(scip, scenariotree, scenarios[stageindex][i]) );

      /* building the tree below the recently added child */
      if( stage < numstages - 1 )
      {
         STOSCENARIO* child = getScenarioChild((*scenariotree), getScenarioNChildren((*scenariotree)) - 1);
         SCIP_CALL( buildScenarioTree(scip, &child, scenarios, numscenarios, numstages, stage + 1) );
      }
   }

   return SCIP_OKAY;
}


/* adds the scenarios to the reader data */
static
SCIP_RETCODE addScenariosToReaderdata(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< the reader data */
   STOSCENARIO***        scenarios,          /**< the array of scenarios */
   int*                  numscenarios,       /**< the number of scenarios per stage */
   int                   numscenariostages   /**< the number of stages for which scenarios were collected */
   )
{
   int i;

   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(scenarios != NULL);
   assert(numscenariostages == SCIPtimGetNStages(scip) - 1);

   SCIP_CALL( buildScenarioTree(scip, &readerdata->scenariotree, scenarios, numscenarios, numscenariostages, 0) );

   /* setting the number of scenarios per stage in the TIME reader data */
   for( i = 0; i < numscenariostages; i++ )
      readerdata->numscenarios += numscenarios[i];

   return SCIP_OKAY;
}


/** finds a scenario with a given name */
static
STOSCENARIO* findScenarioInTree(
   STOSCENARIO*          scenariotree,       /**< the scenario tree to search */
   const char*           scenname            /**< the name of the scenario to search */
   )
{
   STOSCENARIO* retscen;
   int i;

   if( strcmp(getScenarioName(scenariotree), scenname) == 0 )
      return scenariotree;
   else
   {
      retscen = NULL;
      for( i = 0; i < getScenarioNChildren(scenariotree); i++ )
      {
         retscen = findScenarioInTree(scenariotree->children[i], scenname);
         if( retscen != NULL )
            return retscen;
      }
   }

   return NULL;
}


/** inserts a scenario into the reader data scenario tree */
static
SCIP_RETCODE insertScenarioInReaderdata(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< the reader data */
   STOSCENARIO*          scenario,           /**< the scenario to insert in the scenario tree */
   char*                 parentname          /**< the parent scenario for the inserting scenario */
   )
{
   STOSCENARIO* parentscen;

   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(scenario != NULL);

   /* searching for the parent scenario in the tree */
   parentscen = findScenarioInTree(readerdata->scenariotree, parentname);

   /* adding the scenario as a child of the parent scenario */
   SCIP_CALL( scenarioAddChild(scip, &parentscen, scenario) );

   readerdata->numscenarios++;

   return SCIP_OKAY;
}


/* builds the scenarios from the blocks for a given stage */
static
SCIP_RETCODE buildScenariosFromBlocks(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO***        blocks,             /**< the block that form the scenarios */
   STOSCENARIO***        scenarios,          /**< the array to store the scenarios */
   STOSCENARIO***        blocksforscen,      /**< the blocks that will form the scenario */
   int*                  numblocksforscen,   /**< the number of blocks that form the scenario */
   int                   numblocks,          /**< the number of blocks */
   int*                  numblocksperblock,  /**< the number of blocks for a given block */
   int*                  numscenarios,       /**< the number of scenarios */
   int*                  scenariossize,      /**< the size of scenarios array */
   const char*           stage,              /**< the stage for this scenario */
   int                   stagenum,           /**< the number of the stage */
   int                   blocknum            /**< the block number */
   )
{
   SCIP_Bool processed;
   int i;
   int j;

   assert(scip != NULL);
   assert(blocks != NULL);
   assert(scenarios != NULL);
   assert(blocksforscen != NULL);

   processed = FALSE;
   i = blocknum + 1;
   while( !processed && i < numblocks )
   {
      /* it is only necessary to process the next block in the list the belongs to the given stage. */
      if( strcmp(getScenarioStageName(scip, blocks[i][0]), stage) == 0 )
      {
         processed = TRUE;

         for( j = 0; j < numblocksperblock[i]; j++ )
         {
            /* adding the blocks that will build the scenario */
            (*blocksforscen)[(*numblocksforscen)] = blocks[i][j];
            (*numblocksforscen)++;
            SCIP_CALL( buildScenariosFromBlocks(scip, blocks, scenarios, blocksforscen, numblocksforscen, numblocks,
                  numblocksperblock, numscenarios, scenariossize, stage, stagenum + 1, i)  );

            /* the last block needs to be removed so that a new block can be used in its place */
            (*numblocksforscen)--;
         }
      }
      else
      {
         /* the index is only incremented if no block is processed. This is necessary because the value of i is used in
          * the next if statement for identifying whether all blocks have been processed.
          */
         i++;
      }
   }

   /* when all blocks have been inspected, then it is possible to build the scenario */
   if( i == numblocks )
   {
      char scenarioname[SCIP_MAXSTRLEN];

      /* ensuring the correct amount of memory is available */
      if( (*numscenarios) + 1 > (*scenariossize) )
      {
         int newsize;
         newsize = SCIPcalcMemGrowSize(scip, (*numscenarios) + 1);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, scenarios, (*scenariossize), newsize) );
         (*scenariossize) = newsize;
      }

      SCIP_CALL( createScenarioData(scip, &(*scenarios)[(*numscenarios)]) );

      /* setting the scenario name */
      (void) SCIPsnprintf(scenarioname, SCIP_MAXSTRLEN, "Scenario_%s_%d", stage, (*numscenarios));
      SCIP_CALL( setScenarioName(scip, (*scenarios)[(*numscenarios)], scenarioname) );
      SCIP_CALL( setScenarioStageName(scip, (*scenarios)[(*numscenarios)], stage) );
      SCIP_CALL( setScenarioNum(scip, (*scenarios)[(*numscenarios)], (*numscenarios)) );
      SCIP_CALL( setScenarioStageNum(scip, (*scenarios)[(*numscenarios)], stagenum) );

      /* if there is only a single block for the scenario, then we simply copy the block.
       * Otherwise, the blocks are merged into a single scenario */
      if( (*numblocksforscen) == 1 )
         SCIP_CALL( copyScenario(scip, (*blocksforscen)[0], &(*scenarios)[(*numscenarios)], FALSE) );
      else
      {
         SCIP_CALL( copyScenario(scip, (*blocksforscen)[0], &(*scenarios)[(*numscenarios)], FALSE) );
         for( i = 1; i < (*numblocksforscen); i++ )
            SCIP_CALL( mergeScenarios(scip, (*blocksforscen)[i], &(*scenarios)[(*numscenarios)]) );
      }

      (*numscenarios)++;
   }

   return SCIP_OKAY;
}


/* creates the scenarios from the blocks */
static
SCIP_RETCODE createScenariosFromBlocks(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< the reader data */
   STOSCENARIO***        blocks,             /**< the block that form the scenarios */
   int                   numblocks,          /**< the number of blocks */
   int*                  numblocksperblock,  /**< the number of blocks for each block type */
   int                   numstages           /**< the number of stages */
   )
{
   STOSCENARIO*** scenarios;
   STOSCENARIO** blocksforscen;
   int* numscenarios;
   int* scenariossize;
   int numblocksforscen;
   int stagenum;
   char periods[SCIP_MAXSTRLEN];
   int i;
   int j;

   assert(scip != NULL);
   assert(blocks != NULL);

   /* allocating the memory for the scenarios array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scenarios, numstages) );
   SCIP_CALL( SCIPallocBufferArray(scip, &numscenarios, numstages) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scenariossize, numstages) );
   for( i = 0; i < numstages; i++ )
   {
      scenariossize[i] = STO_DEFAULT_BLOCKARRAYSIZE;
      numscenarios[i] = 0;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scenarios[i], scenariossize[i]) );
   }

   /* allocating the memory for the block for scenario array */
   SCIP_CALL( SCIPallocBufferArray(scip, &blocksforscen, numblocks) );

   (void) SCIPsnprintf(periods, SCIP_MAXSTRLEN, "");

   stagenum = 0;
   for( i = 0; i < numblocks; i++ )
   {
      numblocksforscen = 0;
      if( strstr(periods, getScenarioStageName(scip, blocks[i][0])) == NULL )
      {
         /* recording the stage name as processed */
         (void) SCIPsnprintf(periods, SCIP_MAXSTRLEN, "%s_%s", periods, getScenarioStageName(scip, blocks[i][0]));

         SCIP_CALL( buildScenariosFromBlocks(scip, blocks, &scenarios[stagenum], &blocksforscen, &numblocksforscen,
               numblocks, numblocksperblock, &numscenarios[stagenum], &scenariossize[stagenum],
               getScenarioStageName(scip, blocks[i][0]), stagenum, i - 1) );

         stagenum++;
      }
   }

   /* adding the scenarios to the reader data */
   SCIP_CALL( setScenarioNum(scip, readerdata->scenariotree, 0) );
   SCIP_CALL( setScenarioStageNum(scip, readerdata->scenariotree, 0) );
   SCIP_CALL( addScenariosToReaderdata(scip, readerdata, scenarios, numscenarios, numstages) );

   SCIPfreeBufferArray(scip, &blocksforscen);
   for( i = numstages - 1; i >= 0; i-- )
   {
      for( j = numscenarios[i] - 1; j >= 0; j-- )
         SCIP_CALL( freeScenarioTree(scip, &scenarios[i][j]) );
      SCIPfreeBlockMemoryArray(scip, &scenarios[i], scenariossize[i]);
   }
   SCIPfreeBufferArray(scip, &scenariossize);
   SCIPfreeBufferArray(scip, &numscenarios);
   SCIPfreeBlockMemoryArray(scip, &scenarios, numstages);

   return SCIP_OKAY;
}

/** creates the reader data */
static
SCIP_RETCODE createReaderdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);

   /* creating the initial scenario */
   SCIP_CALL( createScenarioData(scip, &readerdata->scenariotree) );

   /* setting the scenario name and stage name */
   SCIP_CALL( setScenarioName(scip, readerdata->scenariotree, "ROOT") );
   SCIP_CALL( setScenarioStageName(scip, readerdata->scenariotree, SCIPtimGetStageName(scip, 0)) );

   return SCIP_OKAY;
}

/** frees the reader data */
static
SCIP_RETCODE freeReaderdata(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);

   /* freeing the scenario tree */
   if( readerdata->scenariotree != NULL )
      SCIP_CALL( freeScenarioTree(scip, &readerdata->scenariotree) );

   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** creates the sto input structure */
static
SCIP_RETCODE stoinputCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   STOINPUT**            stoi,               /**< sto input structure */
   SCIP_FILE*            fp                  /**< file object for the input file */
   )
{
   assert(stoi != NULL);
   assert(fp != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, stoi) );

   (*stoi)->section     = STO_STOCH;
   (*stoi)->stochinfotype = STO_STOCHINFO_NONE;
   (*stoi)->fp          = fp;
   (*stoi)->lineno      = 0;
   (*stoi)->haserror    = FALSE;
   (*stoi)->buf     [0] = '\0';
   (*stoi)->probname[0] = '\0';
   (*stoi)->stochtype[0] = '\0';
   (*stoi)->f0          = NULL;
   (*stoi)->f1          = NULL;
   (*stoi)->f2          = NULL;
   (*stoi)->f3          = NULL;
   (*stoi)->f4          = NULL;
   (*stoi)->f5          = NULL;
   (*stoi)->f6          = NULL;

   return SCIP_OKAY;
}

/** free the sto input structure */
static
void stoinputFree(
   SCIP*                 scip,               /**< SCIP data structure */
   STOINPUT**            stoi                /**< sto input structure */
   )
{
   SCIPfreeBlockMemory(scip, stoi);
}

/** returns the current section */
static
STOSECTION stoinputSection(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->section;
}

/** returns the stochastic information type */
static
STOSTOCHINFO stoinputStochInfoType(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->stochinfotype;
}

/** return the current value of field 0 */
static
const char* stoinputField0(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->f0;
}

/** return the current value of field 1 */
static
const char* stoinputField1(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->f1;
}

/** return the current value of field 2 */
static
const char* stoinputField2(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->f2;
}

/** return the current value of field 3 */
static
const char* stoinputField3(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->f3;
}

/** return the current value of field 4 */
static
const char* stoinputField4(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->f4;
}

/** return the current value of field 5 */
static
const char* stoinputField5(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->f5;
}

/** return the current value of field 6 */
static
const char* stoinputField6(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->f6;
}

/** returns if an error was detected */
static
SCIP_Bool stoinputHasError(
   const STOINPUT*       stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   return stoi->haserror;
}

/** set the section in the sto input structure to given section */
static
void stoinputSetSection(
   STOINPUT*             stoi,               /**< sto input structure */
   STOSECTION            section             /**< section that is set */
   )
{
   assert(stoi != NULL);

   stoi->section = section;
}

/** set the stochastic info type in the sto input structure */
static
void stoinputSetStochInfoType(
   STOINPUT*             stoi,               /**< sto input structure */
   STOSTOCHINFO          stochinfotype       /**< the stochastic infomation type */
   )
{
   assert(stoi != NULL);

   stoi->stochinfotype = stochinfotype;
}

/** set the problem name in the sto input structure to given problem name */
static
void stoinputSetProbname(
   STOINPUT*             stoi,               /**< sto input structure */
   const char*           probname            /**< name of the problem to set */
   )
{
   assert(stoi     != NULL);
   assert(probname != NULL);
   assert(strlen(probname) < sizeof(stoi->probname));

   (void)SCIPmemccpy(stoi->probname, probname, '\0', STO_MAX_NAMELEN - 1);
}

/** set the type name in the sto input structure to given objective name */
static
void stoinputSetStochtype(
   STOINPUT*             stoi,               /**< sto input structure */
   const char*           stochtype           /**< name of the scenario type */
   )
{
   assert(stoi != NULL);
   assert(stochtype != NULL);
   assert(strlen(stochtype) < sizeof(stoi->stochtype));

   (void)SCIPmemccpy(stoi->stochtype, stochtype, '\0', STO_MAX_NAMELEN - 1);
}

static
void stoinputSyntaxerror(
   STOINPUT*             stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   SCIPerrorMessage("Syntax error in line %d\n", stoi->lineno);
   stoi->section  = STO_ENDATA;
   stoi->haserror = TRUE;
}

/** fill the line from \p pos up to column 80 with blanks. */
static
void clearFrom(
   char*                 buf,                /**< buffer to clear */
   unsigned int          pos                 /**< position to start the clearing process */
   )
{
   unsigned int i;

   for(i = pos; i < 80; i++)
      buf[i] = BLANK;
   buf[80] = '\0';
}

/** read a sto format data line and parse the fields. */
static
SCIP_Bool stoinputReadLine(
   STOINPUT*             stoi                /**< sto input structure */
   )
{
   unsigned int len;
   unsigned int i;
   char* s;
   SCIP_Bool is_marker;
   SCIP_Bool is_empty;
   char* nexttok;

   do
   {
      stoi->f0 = stoi->f1 = stoi->f2 = stoi->f3 = stoi->f4 = stoi->f5 = stoi->f6 = 0;
      is_marker = FALSE;

      /* Read until we have not a comment line. */
      do
      {
         stoi->buf[STO_MAX_LINELEN-1] = '\0';
         if( NULL == SCIPfgets(stoi->buf, (int) sizeof(stoi->buf), stoi->fp) )
            return FALSE;
         stoi->lineno++;
      }
      while( *stoi->buf == '*' );   /* coverity[a_loop_bound] */

      /* Normalize line */
      len = (unsigned int) strlen(stoi->buf);

      for( i = 0; i < len; i++ )
      {
         if( (stoi->buf[i] == '\t') || (stoi->buf[i] == '\n') || (stoi->buf[i] == '\r') )
            stoi->buf[i] = BLANK;
      }

      if( len < 80 )
         clearFrom(stoi->buf, len);

      SCIPdebugMessage("line %d: <%s>\n", stoi->lineno, stoi->buf);

      assert(strlen(stoi->buf) >= 80);

      /* Look for new section */
      if( *stoi->buf != BLANK )
      {
         stoi->f0 = SCIPstrtok(&stoi->buf[0], " ", &nexttok);

         assert(stoi->f0 != 0);

         stoi->f1 = SCIPstrtok(NULL, " ", &nexttok);

         return TRUE;
      }

      s = &stoi->buf[1];

      /* At this point it is not clear if we have a indicator field.
       * If there is none (e.g. empty) f1 will be the first name field.
       * If there is one, f2 will be the first name field.
       *
       * Initially comment marks '$' are only allowed in the beginning
       * of the 2nd and 3rd name field. We test all fields but the first.
       * This makes no difference, since if the $ is at the start of a value
       * field, the line will be erroneous anyway.
       */
      do
      {
         if( NULL == (stoi->f1 = SCIPstrtok(s, " ", &nexttok)) )
            break;

         if( (NULL == (stoi->f2 = SCIPstrtok(NULL, " ", &nexttok))) || (*stoi->f2 == '$') )
         {
            stoi->f2 = 0;
            break;
         }

         if( (NULL == (stoi->f3 = SCIPstrtok(NULL, " ", &nexttok))) || (*stoi->f3 == '$') )
         {
            stoi->f3 = 0;
            break;
         }

         if( (NULL == (stoi->f4 = SCIPstrtok(NULL, " ", &nexttok))) || (*stoi->f4 == '$') )
         {
            stoi->f4 = 0;
            break;
         }

         if( (NULL == (stoi->f5 = SCIPstrtok(NULL, " ", &nexttok))) || (*stoi->f5 == '$') )
         {
            stoi->f5 = 0;
            break;
         }

         if( (NULL == (stoi->f6 = SCIPstrtok(NULL, " ", &nexttok))) || (*stoi->f6 == '$') )
            stoi->f6 = 0;
      }
      while( FALSE );

      /* check for empty lines */
      is_empty = (stoi->f0 == NULL && stoi->f1 == NULL);
   }
   while( is_marker || is_empty );

   return TRUE;
}

/** Process STOCH section. */
static
SCIP_RETCODE readStoch(
   SCIP*                 scip,               /**< SCIP data structure */
   STOINPUT*             stoi                /**< sto input structure */
   )
{
   assert(stoi != NULL);

   SCIPdebugMsg(scip, "read problem name\n");

   /* This has to be the Line with the NAME section. */
   if( !stoinputReadLine(stoi) || stoinputField0(stoi) == NULL || strcmp(stoinputField0(stoi), "STOCH") )
   {
      stoinputSyntaxerror(stoi);
      return SCIP_OKAY;
   }

   /* Sometimes the name is omitted. */
   stoinputSetProbname(stoi, (stoinputField1(stoi) == 0) ? "_STO_" : stoinputField1(stoi));

   /* This hat to be a new section */
   /* coverity[tainted_data] */
   if( !stoinputReadLine(stoi) || (stoinputField0(stoi) == NULL) )
   {
      stoinputSyntaxerror(stoi);
      return SCIP_OKAY;
   }

   /* setting the stochatic information section */
   if( !strncmp(stoinputField0(stoi), "BLOCKS", 6) )
      stoinputSetSection(stoi, STO_BLOCKS);
   else if( !strncmp(stoinputField0(stoi), "SCENARIOS", 9) )
      stoinputSetSection(stoi, STO_SCENARIOS);
   else if( !strncmp(stoinputField0(stoi), "INDEP", 5) )
      stoinputSetSection(stoi, STO_INDEP);
   else
   {
      stoinputSyntaxerror(stoi);
      return SCIP_OKAY;
   }

   /* setting the stochastic information type */
   if( !strncmp(stoinputField1(stoi), "DISCRETE", 8) )
      stoinputSetStochInfoType(stoi, STO_STOCHINFO_DISCRETE);
   else if( !strncmp(stoinputField1(stoi), "UNIFORM", 7) )
      stoinputSetStochInfoType(stoi, STO_STOCHINFO_UNIFORM);
   else if( !strncmp(stoinputField1(stoi), "NORMAL", 6) )
      stoinputSetStochInfoType(stoi, STO_STOCHINFO_NORMAL);
   else if( !strncmp(stoinputField1(stoi), "SUB", 3) )
      stoinputSetStochInfoType(stoi, STO_STOCHINFO_SUB);
   else if( !strncmp(stoinputField1(stoi), "LINTR", 5) )
      stoinputSetStochInfoType(stoi, STO_STOCHINFO_LINTR);
   else
   {
      stoinputSyntaxerror(stoi);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** Process BLOCKS section. */
static
SCIP_RETCODE readBlocks(
   STOINPUT*             stoi,               /**< sto input structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   STOSCENARIO*** blocks;
   int numblocks;
   int* numblocksperblock;
   int blockssize;
   int* blocksperblocksize;
   char BL[] = "BL";
   int blocknum;
   int blockindex;
   int i;
   int j;
   char stagenames[SCIP_MAXSTRLEN];
   int numstages;

   SCIPdebugMsg(scip, "read Blocks\n");

   /* This has to be the Line with the name. */
   if( stoinputField1(stoi) == NULL )
   {
      stoinputSyntaxerror(stoi);
      return SCIP_OKAY;
   }

   stoinputSetStochtype(stoi, stoinputField1(stoi));

   /* initializing the block data */
   numblocks = 0;
   blockssize = STO_DEFAULT_ARRAYSIZE;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blocks, STO_DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &numblocksperblock, STO_DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blocksperblocksize, STO_DEFAULT_ARRAYSIZE) );

   blockindex = 0;
   blocknum = 0;

   /* initializing the stage names record */
   numstages = 0;
   (void) SCIPsnprintf(stagenames, SCIP_MAXSTRLEN, "");

   /* coverity[tainted_data] */
   while( stoinputReadLine(stoi) )
   {
      if( stoinputField0(stoi) != NULL )
      {
         if( !strcmp(stoinputField0(stoi), "BLOCKS") )
         {
            stoinputSetSection(stoi, STO_BLOCKS);
            if( strcmp(stoinputField1(stoi), "DISCRETE") )
            {
               SCIPerrorMessage("Sorry, %s blocks stucture is not currently supported.\n", stoinputField1(stoi));
               SCIPerrorMessage("Only DISCRETE blocks are supported.\n");
               goto TERMINATE;
            }
         }
         else if( !strcmp(stoinputField0(stoi), "ENDATA") )
         {
            SCIP_CALL( createScenariosFromBlocks(scip, readerdata, blocks, numblocks, numblocksperblock, numstages) );
            stoinputSetSection(stoi, STO_ENDATA);
         }
         else
            stoinputSyntaxerror(stoi);

         goto TERMINATE;
      }

      if( strcmp(stoinputField1(stoi), BL) == 0 )
      {
         SCIP_Bool foundblock = FALSE;

         /* checking whether the stage has been added previously */
         if( strstr(stagenames, stoinputField3(stoi)) == NULL )
         {
            /* recording the stage name as processed */
            (void) SCIPsnprintf(stagenames, SCIP_MAXSTRLEN, "%s_%s", stagenames, stoinputField3(stoi));
            numstages++;
         }

         /* determining whether a block name has previously been added */
         for( i = 0; i < numblocks; i++ )
         {
            if( strcmp(getScenarioName(blocks[i][0]), stoinputField2(stoi)) == 0 )
            {
               foundblock = TRUE;
               break;
            }
         }
         blocknum = i;

         /* if the block is found, then the memory for the blocks array must be ensured */
         if( foundblock )
         {
            /* ensuring enough memory is available for the blocks */
            if( numblocksperblock[blocknum] + 1 > blocksperblocksize[blocknum] )
            {
               int newsize;
               newsize = SCIPcalcMemGrowSize(scip, numblocksperblock[blocknum] + 1);
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &blocks[blocknum], blocksperblocksize[blocknum], newsize) ); /*lint !e866*/
               blocksperblocksize[blocknum] = newsize;
            }
         }
         else
         {
            /* ensuring enough memory is available for the blocks */
            if( numblocks + 1 > blockssize )
            {
               int newsize;
               newsize = SCIPcalcMemGrowSize(scip, numblocks + 1);
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &blocks, blockssize, newsize) );
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &numblocksperblock, blockssize, newsize) );
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &blocksperblocksize, blockssize, newsize) );
               blockssize = newsize;
            }

            blocksperblocksize[blocknum] = STO_DEFAULT_BLOCKARRAYSIZE;
            numblocksperblock[blocknum] = 0;
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blocks[blocknum], blocksperblocksize[blocknum]) );
         }

         blockindex = numblocksperblock[blocknum];

         /* creating the scenario data structure */
         SCIP_CALL( createScenarioData(scip, &blocks[blocknum][blockindex]) );

         SCIP_CALL( setScenarioName(scip, blocks[blocknum][blockindex], stoinputField2(stoi)) );
         SCIP_CALL( setScenarioStageName(scip, blocks[blocknum][blockindex], stoinputField3(stoi)) );
         SCIP_CALL( setScenarioProbability(scip, blocks[blocknum][blockindex], atof(stoinputField4(stoi))) );
         numblocksperblock[blocknum]++;

         if( !foundblock )
            numblocks++;
      }
      else
      {
         SCIP_CALL( addScenarioEntry(scip, blocks[blocknum][blockindex], stoinputField2(stoi), stoinputField1(stoi),
               atof(stoinputField3(stoi))) );
      }
   }
   stoinputSyntaxerror(stoi);

TERMINATE:

   /* releasing the scenario data */
   for( i = numblocks - 1; i >= 0; i-- )
   {
      for( j = numblocksperblock[i] - 1; j >= 0; j-- )
         SCIP_CALL( freeScenarioTree(scip, &blocks[i][j]) );
   }

   for( i = numblocks - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &blocks[i], blocksperblocksize[i]);
   SCIPfreeBlockMemoryArray(scip, &blocksperblocksize, blockssize);
   SCIPfreeBlockMemoryArray(scip, &numblocksperblock, blockssize);
   SCIPfreeBlockMemoryArray(scip, &blocks, blockssize);

   return SCIP_OKAY;
}


/** Process SCENARIOS section. */
static
SCIP_RETCODE readScenarios(
   STOINPUT*             stoi,               /**< sto input structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   STOSCENARIO* scenario;
   char SC[] = "SC";
   char wrongroot[] = "\'ROOT\'";
   char parentname[SCIP_MAXSTRLEN];
   char scennames[SCIP_MAXSTRLEN];
   char tmpname[SCIP_MAXSTRLEN];
   int numscenarios;
   SCIP_Bool addscenario;

   SCIPdebugMsg(scip, "read SCENARIOS\n");

   /* This has to be the Line with the name. */
   if( stoinputField1(stoi) == NULL )
   {
      stoinputSyntaxerror(stoi);
      return SCIP_OKAY;
   }

   stoinputSetStochtype(stoi, stoinputField1(stoi));

   /* initializing the scen names record */
   numscenarios = 0;
   (void) SCIPsnprintf(scennames, SCIP_MAXSTRLEN, "ROOT");

   scenario = NULL;
   addscenario = FALSE;

   /* initializing the root scenario in the reader data */
   SCIP_CALL( setScenarioNum(scip, readerdata->scenariotree, 0) );
   SCIP_CALL( setScenarioStageNum(scip, readerdata->scenariotree, 0) );

   /* coverity[tainted_data] */
   while( stoinputReadLine(stoi) )
   {
      if( stoinputField0(stoi) != NULL )
      {
         /* if a scenario has been created that needs to be added to the scenario tree */
         if( addscenario )
         {
            SCIP_CALL( insertScenarioInReaderdata(scip, readerdata, scenario, parentname) );

            /* freeing the scenario */
            SCIP_CALL( freeScenarioTree(scip, &scenario) );
         }

         if( !strcmp(stoinputField0(stoi), "SCENARIOS") )
         {
            stoinputSetSection(stoi, STO_SCENARIOS);
            if( strcmp(stoinputField1(stoi), "DISCRETE") )
            {
               SCIPerrorMessage("Sorry, %s scenarios is not currently supported.\n", stoinputField1(stoi));
               SCIPerrorMessage("Only DISCRETE scenarios are supported.\n");
               goto TERMINATE;
            }
         }
         else if( !strcmp(stoinputField0(stoi), "ENDATA") )
            stoinputSetSection(stoi, STO_ENDATA);
         else
            stoinputSyntaxerror(stoi);

         goto TERMINATE;
      }

      if( strcmp(stoinputField1(stoi), SC) == 0 )
      {
         int stagenum;

         /* if a scenario has been created that needs to be added to the scenario tree */
         if( addscenario )
         {
            SCIP_CALL( insertScenarioInReaderdata(scip, readerdata, scenario, parentname) );

            /* freeing the scenario */
            SCIP_CALL( freeScenarioTree(scip, &scenario) );
            assert(scenario == NULL);
         }

         if( strcmp(wrongroot, stoinputField3(stoi)) == 0 )
            (void) SCIPsnprintf(parentname, SCIP_MAXSTRLEN, "%s", "ROOT");
         else
            (void) SCIPsnprintf(parentname, SCIP_MAXSTRLEN, "%s", stoinputField3(stoi));

         /* checking whether the stage has been added previously */
         if( strstr(scennames, stoinputField2(stoi)) == NULL )
         {
            /* recording the stage name as processed */
            (void) SCIPsnprintf(tmpname, SCIP_MAXSTRLEN, "%s_%s", scennames, stoinputField2(stoi));
            (void) SCIPsnprintf(scennames, SCIP_MAXSTRLEN, "%s", tmpname);
         }

         /* checking whether the "common" scenario has been added yet */
         if( strstr(scennames, parentname) == NULL )
         {
            SCIPerrorMessage("Scenario <%s> needs to be read before scenario <%s>\n", parentname, stoinputField2(stoi));
            stoinputSyntaxerror(stoi);
            goto TERMINATE;
         }

         /* the "common" scenario has been added before, so a child can be added to the scenario tree */
         SCIP_CALL( createScenarioData(scip, &scenario) );

         SCIP_CALL( setScenarioName(scip, scenario, stoinputField2(stoi)) );
         SCIP_CALL( setScenarioStageName(scip, scenario, stoinputField5(stoi)) );
         SCIP_CALL( setScenarioNum(scip, scenario, numscenarios) );

         stagenum = SCIPtimFindStage(scip, stoinputField5(stoi));
         if( stagenum < 0 )
         {
            stoinputSyntaxerror(stoi);
            goto TERMINATE;
         }
         SCIP_CALL( setScenarioStageNum(scip, scenario, stagenum) );
         SCIP_CALL( setScenarioProbability(scip, scenario, atof(stoinputField4(stoi))) );
         if( stoinputField6(stoi) != NULL )
         {
            SCIP_CALL( setScenarioLowerbound(scip, scenario, atof(stoinputField6(stoi))) );
         }

         numscenarios++;
         addscenario = TRUE;
      }
      else if( addscenario )
      {
         SCIP_CALL( addScenarioEntry(scip, scenario, stoinputField2(stoi), stoinputField1(stoi),
               atof(stoinputField3(stoi))) );
      }
   }
   stoinputSyntaxerror(stoi);

TERMINATE:

   return SCIP_OKAY;
}


/** Process INDEP section. */
static
SCIP_RETCODE readIndep(
   STOINPUT*             stoi,               /**< sto input structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   STOSCENARIO*** blocks;
   int numblocks;
   int* numblocksperblock;
   int blockssize;
   int* blocksperblocksize;
   int blocknum;
   int blockindex;
   int i;
   int j;
   char stagenames[SCIP_MAXSTRLEN];
   int numstages;
   SCIP_Bool foundblock;

   SCIP_Real probability;
   char currstagename[SCIP_MAXSTRLEN];

   SCIPdebugMsg(scip, "read Indep\n");

   /* This has to be the Line with the name. */
   if( stoinputField1(stoi) == NULL )
   {
      stoinputSyntaxerror(stoi);
      return SCIP_OKAY;
   }

   stoinputSetStochtype(stoi, stoinputField1(stoi));

   /* initializing the block data */
   numblocks = 0;
   blockssize = STO_DEFAULT_ARRAYSIZE;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blocks, STO_DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &numblocksperblock, STO_DEFAULT_ARRAYSIZE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blocksperblocksize, STO_DEFAULT_ARRAYSIZE) );

   /* initializing the stage names record */
   numstages = 0;
   (void) SCIPsnprintf(stagenames, SCIP_MAXSTRLEN, "");

   while( stoinputReadLine(stoi) )
   {
      if( stoinputField0(stoi) != NULL )
      {
         if( !strcmp(stoinputField0(stoi), "INDEP") )
         {
            stoinputSetSection(stoi, STO_INDEP);
         }
         else if( !strcmp(stoinputField0(stoi), "ENDATA") )
         {
            SCIP_CALL( createScenariosFromBlocks(scip, readerdata, blocks, numblocks, numblocksperblock, numstages) );
            stoinputSetSection(stoi, STO_ENDATA);
         }
         else
            stoinputSyntaxerror(stoi);

         goto TERMINATE;
      }

      /* if the 5th input is NULL, then the 4th input is the probability. Otherwise, the 4th input is the stage name and
       * the 5th input is the probability. The stage name is redundant information, but sometimes included for more
       * information.
       */
      if( stoinputField5(stoi) == NULL )
      {
         probability = atof(stoinputField4(stoi));
         (void) SCIPsnprintf(currstagename, SCIP_MAXSTRLEN, "%s", SCIPtimConsGetStageName(scip, stoinputField2(stoi)));
      }
      else
      {
         probability = atof(stoinputField5(stoi));
         (void) SCIPsnprintf(currstagename, SCIP_MAXSTRLEN, "%s", stoinputField4(stoi));
      }

      /* checking whether the stage has been added previously */
      if( strstr(stagenames, currstagename) == NULL )
      {
         /* recording the stage name as processed */
         (void) SCIPsnprintf(stagenames, SCIP_MAXSTRLEN, "%s_%s", stagenames, currstagename);

         numstages++;
      }

      foundblock = FALSE;

      /* determining whether a block name has previously been added */
      for( i = 0; i < numblocks; i++ )
      {
         if( strcmp(getScenarioName(blocks[i][0]), stoinputField2(stoi)) == 0 )
         {
            foundblock = TRUE;
            break;
         }
      }
      blocknum = i;

      /* if the block is found, then the memory for the blocks array must be ensured */
      if( foundblock )
      {
         /* ensuring enough memory is available for the blocks */
         if( numblocksperblock[blocknum] + 1 > blocksperblocksize[blocknum] )
         {
            int newsize;
            newsize = SCIPcalcMemGrowSize(scip, numblocksperblock[blocknum] + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &blocks[blocknum], blocksperblocksize[blocknum], newsize) );   /*lint !e866*/
            blocksperblocksize[blocknum] = newsize;
         }
      }
      else
      {
         /* ensuring enough memory is available for the blocks */
         if( numblocks + 1 > blockssize )
         {
            int newsize;
            newsize = SCIPcalcMemGrowSize(scip, numblocks + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &blocks, blockssize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &numblocksperblock, blockssize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &blocksperblocksize, blockssize, newsize) );
            blockssize = newsize;
         }

         blocksperblocksize[blocknum] = STO_DEFAULT_BLOCKARRAYSIZE;
         numblocksperblock[blocknum] = 0;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blocks[blocknum], blocksperblocksize[blocknum]) );
      }

      blockindex = numblocksperblock[blocknum];

      /* creating the scenario data structure */
      SCIP_CALL( createScenarioData(scip, &blocks[blocknum][blockindex]) );

      SCIP_CALL( setScenarioName(scip, blocks[blocknum][blockindex], stoinputField2(stoi)) );
      SCIP_CALL( setScenarioStageName(scip, blocks[blocknum][blockindex], currstagename) );
      SCIP_CALL( setScenarioProbability(scip, blocks[blocknum][blockindex], probability) );
      numblocksperblock[blocknum]++;

      if( !foundblock )
         numblocks++;

      SCIP_CALL( addScenarioEntry(scip, blocks[blocknum][blockindex], stoinputField2(stoi), stoinputField1(stoi),
            atof(stoinputField3(stoi))) );
   }
   stoinputSyntaxerror(stoi);

TERMINATE:

   /* releasing the scenario data */
   for( i = numblocks - 1; i >= 0; i-- )
   {
      for( j = numblocksperblock[i] - 1; j >= 0; j-- )
         SCIP_CALL( freeScenarioTree(scip, &blocks[i][j]) );
   }

   for( i = numblocks - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &blocks[i], blocksperblocksize[i]);
   SCIPfreeBlockMemoryArray(scip, &blocksperblocksize, blockssize);
   SCIPfreeBlockMemoryArray(scip, &numblocksperblock, blockssize);
   SCIPfreeBlockMemoryArray(scip, &blocks, blockssize);

   return SCIP_OKAY;
}


/** computes the probability of a scenario */
static
SCIP_Real computeScenarioProbability(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario            /**< the current scenario */
   )
{
   STOSCENARIO* checkscen;
   SCIP_Real probability;

   assert(scip != NULL);
   assert(scenario != NULL);

   /* computing the probability for the scenario */
   checkscen = scenario;
   probability = 1;
   while( checkscen != NULL )
   {
      probability *= getScenarioProbability(scip, checkscen);
      checkscen = getScenarioParent(checkscen);
   }

   return probability;
}

/** gets the variable name */
static
void getScenarioEntityName(
   char*                 name,               /**< the name to be returned */
   const char*           varname,            /**< the root of the variable name */
   int                   stagenum,           /**< the stage number */
   int                   scenarionum         /**< the scenario number */
   )
{
   if( stagenum < 0 )
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_00_%d", varname, scenarionum);
   else
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d_%d", varname, stagenum, scenarionum);
}


/** add variables to the scenario  */
static
SCIP_RETCODE addScenarioVarsToProb(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the current scenario */
   SCIP_HASHMAP*         varmap,             /**< the variable map from the original to the subproblem variables */
   SCIP_VAR**            vars,               /**< the variables of the core problem associated with this scenario */
   int                   nvars               /**< the number of variables for this scenario */
   )
{
   SCIP_Real probability;
   int i;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(scenario != NULL);
   assert(vars != NULL);

   /* computing the probability for the scenario */
   probability = computeScenarioProbability(scip, scenario);

   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_VARTYPE vartype;

      SCIPdebugMessage("Original problem variable <%s> is being duplicated for scenario %d\n", SCIPvarGetName(vars[i]),
         getScenarioNum(scip, scenario));

      if( SCIPvarIsDeleted(vars[i]) )
         continue;

      obj = SCIPvarGetObj(vars[i])*probability;

      vartype = SCIPvarGetType(vars[i]);

      /* creating a variable as a copy of the original variable. */
      getScenarioEntityName(name, SCIPvarGetName(vars[i]), getScenarioStageNum(scip, scenario), getScenarioNum(scip, scenario));
      SCIP_CALL( SCIPcreateVar(scip, &var, name, SCIPvarGetLbOriginal(vars[i]), SCIPvarGetUbOriginal(vars[i]),
            obj, vartype, SCIPvarIsInitial(vars[i]), SCIPvarIsRemovable(vars[i]), NULL, NULL, NULL,
            NULL, NULL) );

      SCIPdebugMessage("Adding variable <%s>\n", name);

      SCIP_CALL( SCIPaddVar(scip, var) );

      /* inserting the scenario variable into the hashmap */
      SCIP_CALL( SCIPhashmapInsert(varmap, vars[i], var) );

      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}


/** finds the scenario variable to add to a constraint */
static
SCIP_RETCODE findScenarioVar(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the current scenario */
   SCIP_VAR*             consvar,            /**< the variable in the constraint that is being searched for */
   SCIP_VAR**            scenariovar         /**< pointer to return the variable to be added to the constraint */
   )
{
   STOSCENARIO* checkscen;
   char varname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(scenario != NULL);
   assert(consvar != NULL);
   assert(scenariovar != NULL);

   (*scenariovar) = NULL;

   checkscen = scenario;

   /* NOTE: if the variable does not exist, then we need to search the preceding scenarios. In the case of
    * decomposition, then we only check the preceding scenario. As such, a check count is used to limit the number
    * of scenario checks. */
   while( (*scenariovar) == NULL )
   {
      assert(checkscen != NULL);
      if( getScenarioStageNum(scip, checkscen) == 0 )
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", SCIPvarGetName(consvar));
      else
         getScenarioEntityName(varname, SCIPvarGetName(consvar), getScenarioStageNum(scip, checkscen),
            getScenarioNum(scip, checkscen));

      (*scenariovar) = SCIPfindVar(scip, varname);

      checkscen = getScenarioParent(checkscen);
   }

   if( (*scenariovar) == NULL )
   {
      SCIPerrorMessage("There is no scenario variable could be found.\n");
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}


/** create variable for the decomposed scenario */
static
SCIP_RETCODE getScenarioDecompVar(
   SCIP*                 scip,               /**< the SCIP data structure */
   STOSCENARIO*          scenario,           /**< the current scenario */
   SCIP_VAR*             consvar,            /**< the variable in the constraint that is being searched for */
   SCIP_VAR**            scenariovar,        /**< pointer to return the variable to be added to the constraint */
   SCIP_Bool*            varadded            /**< pointer to indicate whether a variable has been added */
   )
{
   STOSCENARIO* checkscen;
   SCIP_VAR* searchvar;
   int checkcount;
   char varname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(scenario != NULL);
   assert(consvar != NULL);

   (*varadded) = FALSE;

   /* finding the scenario that the consvar belongs to */
   checkscen = scenario;
   searchvar = NULL;
   checkcount = 0;
   while( searchvar == NULL && checkcount < 2 )
   {
      assert(checkscen != NULL);
      if( getScenarioStageNum(scip, checkscen) == 0 )
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", SCIPvarGetName(consvar));
      else
         getScenarioEntityName(varname, SCIPvarGetName(consvar), getScenarioStageNum(scip, checkscen),
            getScenarioNum(scip, checkscen));

      /* first checking whether the variable is included in the scenario */
      searchvar = SCIPfindVar(scip, varname);
      if( searchvar != NULL )
      {
         (*scenariovar) = searchvar;
         return SCIP_OKAY;
      }

      searchvar = SCIPfindVar(getScenarioScip(checkscen), varname);

      checkscen = getScenarioParent(checkscen);
      checkcount++;
   }

   if( searchvar != NULL )
   {
      SCIP_VAR* var;
      /* creating a variable as a copy of the original variable. */
      SCIP_CALL( SCIPcreateVar(scip, &var, varname, SCIPvarGetLbOriginal(searchvar), SCIPvarGetUbOriginal(searchvar),
            0.0, SCIPvarGetType(searchvar), SCIPvarIsInitial(searchvar), SCIPvarIsRemovable(searchvar), NULL, NULL,
            NULL, NULL, NULL) );

      SCIP_CALL( SCIPaddVar(scip, var) );

      (*scenariovar) = var;
      (*varadded) = TRUE;
   }

   return SCIP_OKAY;
}


/** adds the constraint to the scenario problem */
static
SCIP_RETCODE addScenarioConsToProb(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 scenarioscip,       /**< the scenario SCIP data structure */
   STOSCENARIO*          scenario,           /**< the current scenario */
   SCIP_HASHMAP*         varmap,             /**< the variable map from the original to the subproblem variables */
   SCIP_CONS**           conss,              /**< the constraints of the core problem associated with this scenario */
   int                   nconss,             /**< the number of constraints for this scenario */
   SCIP_Bool             decomp              /**< is the problem being decomposed */
   )
{
   int i;
   int j;
   char name[SCIP_MAXSTRLEN];
   SCIP_Bool varadded;

   assert(scip != NULL);
   assert(scenarioscip != NULL);
   assert(scenario != NULL);
   assert(conss != NULL);

   /* Add constraints */
   /* NOTE: It is assumed that the problems only have linear constraints */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars = NULL;
      int nconsvars;
      SCIP_Bool success1 = TRUE;
      SCIP_Bool success2 = TRUE;

      if( SCIPconsIsDeleted(conss[i]) )
         continue;

      /* getting the number of variables in the constraints */
      SCIP_CALL( SCIPgetConsNVars(scip, conss[i], &nconsvars, &success1) );

      if( success1 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
         SCIP_CALL( SCIPgetConsVars(scip, conss[i], consvars, nconsvars, &success2) );

         /* If the get variable callback is not implemented for the constraint, then the success flag will be returned
          * as FALSE. In this case, it is not possible to build the stochastic program, so an error will be returned.
          */
         if( !success2 )
         {
            SCIPfreeBufferArrayNull(scip, consvars);
         }
      }

      if( !success1 || !success2 )
      {
         SCIPerrorMessage("It is not possible to copy constraint <%s>. The stochastic program can not be built.\n",
            SCIPconsGetName(conss[i]));

         return SCIP_READERROR;
      }

      assert(consvars != NULL);
      for( j = 0; j < nconsvars; j++ )
      {
         SCIP_VAR* scenariovar;

         scenariovar = NULL;

         varadded = FALSE;

         if( decomp )
            SCIP_CALL( getScenarioDecompVar(scenarioscip, scenario, consvars[j], &scenariovar, &varadded) );
         else
            SCIP_CALL( findScenarioVar(scenarioscip, scenario, consvars[j], &scenariovar) );

         if( scenariovar != NULL )
         {
            /* checking whether the variable is in the variable hashmap. If it doesn't exist, then it is added to the
             * variable hashmap
             */
            if( !SCIPhashmapExists(varmap, consvars[j]) )
            {
               SCIP_CALL( SCIPhashmapInsert(varmap, consvars[j], scenariovar) );
            }
         }

         if( varadded )
         {
            SCIP_CALL( SCIPreleaseVar(scenarioscip, &scenariovar) );
         }
      }

      /* creating a linear constraint as a copy of the original constraint. */
      getScenarioEntityName(name, SCIPconsGetName(conss[i]), getScenarioStageNum(scip, scenario), getScenarioNum(scip, scenario));

      /* copying the constraint from the original SCIP to the stochastic program */
      SCIP_CALL( SCIPgetConsCopy(scip, scenarioscip, conss[i], &cons, SCIPconsGetHdlr(conss[i]), varmap, NULL, name,
            SCIPconsIsInitial(conss[i]), SCIPconsIsSeparated(conss[i]), SCIPconsIsEnforced(conss[i]),
            SCIPconsIsChecked(conss[i]), SCIPconsIsMarkedPropagate(conss[i]), SCIPconsIsLocal(conss[i]),
            SCIPconsIsModifiable(conss[i]), SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]),
            SCIPconsIsStickingAtNode(conss[i]), TRUE, &success1) );

      /* freeing the cons vars buffer array */
      SCIPfreeBufferArray(scip, &consvars);

      /* if the copy failed, then the scenarios can not be created. */
      if( !success1 )
      {
         SCIPerrorMessage("It is not possible to copy constraint <%s>. The stochastic program can not be built.\n",
            SCIPconsGetName(conss[i]));
         return SCIP_READERROR;
      }

      SCIP_CALL( SCIPaddCons(scenarioscip, cons) );
      SCIP_CALL( SCIPreleaseCons(scenarioscip, &cons) );
   }

   return SCIP_OKAY;
}

/** add variables and constraint to problem */
static
SCIP_RETCODE addScenarioVarsAndConsToProb(
   SCIP*                 scip,               /**< the SCIP data structure of master problem */
   STOSCENARIO*          scenario,           /**< the current scenario */
   SCIP_Bool             decomp              /**< is the problem being decomposed */
   )
{
   SCIP* scenarioscip;
   SCIP_BENDERS* benders;
   SCIP_HASHMAP* varmap;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_Real probability;
   int nconss;
   int nvars;
   int nmastervars;
   int nentries;
   int stagenum;
   int i;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(scenario != NULL);

   stagenum = SCIPtimFindStage(scip, getScenarioStageName(scip, scenario));
   if( stagenum < 0 || stagenum >= SCIPtimGetNStages(scip) )
   {
      SCIPerrorMessage("Unable to build stochastic program - stage <%s> was not found\n",
         getScenarioStageName(scip, scenario));
      return SCIP_READERROR;
   }

   SCIPdebugMessage("Creating scenario at stage <%d>. Scenario: %d Stage: %d\n", stagenum, getScenarioNum(scip, scenario),
      getScenarioStageNum(scip, scenario));

   conss = SCIPtimGetStageConss(scip, stagenum);
   nconss = SCIPtimGetStageNConss(scip, stagenum);
   vars = SCIPtimGetStageVars(scip, stagenum);
   nvars = SCIPtimGetStageNVars(scip, stagenum);

   nmastervars = SCIPgetNVars(scip);

   /* this if 0 will be removed when the stochastic reader is merged with the Benders' branch */
   if( decomp )
   {
      SCIP_CALL( SCIPcreate(&scenarioscip) );

      getScenarioEntityName(name, SCIPgetProbName(scip), getScenarioStageNum(scip, scenario), getScenarioNum(scip, scenario));

      /* creating the problem */
      SCIP_CALL( SCIPcreateProbBasic(scenarioscip, name) );

      /* we explicitly enable the use of a debug solution for this main SCIP instance */
      SCIPenableDebugSol(scenarioscip);

      /* include default SCIP plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(scenarioscip) );

      /* activating the Benders' constraint handler for the scenario stages.
       * TODO: consider whether the two-phase method should be activated by default in the scenario stages.
       */
      SCIP_CALL( SCIPsetBoolParam(scenarioscip, "constraints/benders/active", TRUE) );

      /* allocating memory for the subproblems */
      if( getScenarioNChildren(scenario) > 0 )
         SCIP_CALL( createScenarioSubproblemArray(scip, scenario) );
   }
   else
      scenarioscip = scip;

   /* adding the scenarioscip to the scenario */
   setScenarioScip(scenario, scenarioscip);

   /* creating the variable hashmap to copy the constraints */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scenarioscip), nmastervars) );

   /* adding the variables to the scenario */
   SCIP_CALL( addScenarioVarsToProb(scenarioscip, scenario, varmap, vars, nvars) );

   /* adding the constraints to the scenario */
   SCIP_CALL( addScenarioConsToProb(scip, scenarioscip, scenario, varmap, conss, nconss, decomp) );

   /* destroying the hashmap */
   SCIPhashmapFree(&varmap);

   /* add the variables and constraints of the child scenarios */
   for( i = 0; i < getScenarioNChildren(scenario); i++ )
   {
      /* the master SCIP is always passed to the recursive function. The scenario SCIP instances are generated in the
       * function call. */
      SCIP_CALL( addScenarioVarsAndConsToProb(scip, getScenarioChild(scenario, i), decomp) );
      if( decomp )
         addScenarioSubproblem(scenario, getScenarioScip(getScenarioChild(scenario, i)));
   }

   /* adding the Benders' decomposition */
   if( decomp && getScenarioNChildren(scenario) > 0 )
   {
      SCIP_CALL( SCIPcreateBendersDefault(scenarioscip, getScenarioSubproblemArray(scenario), getScenarioNChildren(scenario)) );

      /* getting the default Benders' decomposition */
      benders = SCIPfindBenders(scenarioscip, "default");

      /* updating the lower bounds for the subproblems */
      for( i = 0; i < getScenarioNChildren(scenario); i++ )
         SCIPbendersUpdateSubproblemLowerbound(benders, i,
            getScenarioLowerbound(scenarioscip, getScenarioChild(scenario, i)));
   }

   /* computing the probability for the scenario */
   probability = computeScenarioProbability(scenarioscip, scenario);

   /* change the constraints for the given scenario */
   nentries = getScenarioNEntries(scenario);
   for( i = 0; i < nentries; i++ )
   {
      SCIP_CONS* cons;
      SCIP_VAR* var;
      char RHS[] = "RHS";
      char rhs[] = "rhs";
      char RIGHT[] = "RIGHT";
      char MINI[] = "MINI";
      char obj[] = "obj";
      char OBJ[] = "OBJ";

      /* finding the constraint associated with the row */
      getScenarioEntityName(name, getScenarioEntryRow(scenario, i), getScenarioStageNum(scenarioscip, scenario),
         getScenarioNum(scenarioscip, scenario));
      cons = SCIPfindCons(scenarioscip, name);

      if( strncmp(getScenarioEntryCol(scenario, i), RHS, 3) == 0 ||
         strncmp(getScenarioEntryCol(scenario, i), rhs, 3) == 0 ||
         strcmp(getScenarioEntryCol(scenario, i), RIGHT) == 0 )
      {
         /* if the constraint is NULL, then it is not possible to make any changes to the scenario */
         if( cons == NULL )
         {
            SCIPerrorMessage("There is no constraint <%s> in the current scenario.\n", name);
            return SCIP_READERROR;
         }

         /* if the constraint is an equality constraint, then the LHS must also be changed */
         if( SCIPgetLhsLinear(scenarioscip, cons) >= SCIPgetRhsLinear(scenarioscip, cons) )
         {
            SCIP_CALL( SCIPchgLhsLinear(scenarioscip, cons, getScenarioEntryValue(scenario, i)) );
            SCIP_CALL( SCIPchgRhsLinear(scenarioscip, cons, getScenarioEntryValue(scenario, i)) );
         }
         else if( SCIPisLT(scenarioscip, SCIPgetRhsLinear(scenarioscip, cons), SCIPinfinity(scenarioscip)) )
            SCIP_CALL( SCIPchgRhsLinear(scenarioscip, cons, getScenarioEntryValue(scenario, i)) );
         else if( SCIPisLT(scenarioscip, SCIPgetLhsLinear(scenarioscip, cons), SCIPinfinity(scenarioscip)) )
            SCIP_CALL( SCIPchgLhsLinear(scenarioscip, cons, getScenarioEntryValue(scenario, i)) );
      }
      else if( strstr(getScenarioEntryRow(scenario, i), MINI) != NULL ||
         strstr(getScenarioEntryRow(scenario, i), obj) != NULL ||
         strstr(getScenarioEntryRow(scenario, i), OBJ) != NULL )
      {
         /* finding the variable associated with the column */
         getScenarioEntityName(name, getScenarioEntryCol(scenario, i), getScenarioStageNum(scenarioscip, scenario),
            getScenarioNum(scenarioscip, scenario));
         var = SCIPfindVar(scenarioscip, name);

         /* changing the coefficient for the variable */
         if( var == NULL )
         {
            SCIPerrorMessage("There is no variable <%s> in the current scenario.\n", name);
            return SCIP_READERROR;
         }
         else
         {
            SCIP_CALL( SCIPchgVarObj(scenarioscip, var, getScenarioEntryValue(scenario, i)*probability) );
         }
      }
      else
      {
         /* if the constraint is NULL, then it is not possible to make any changes to the scenario */
         if( cons == NULL )
         {
            SCIPerrorMessage("There is no constraint <%s> in the current scenario.\n", name);
            return SCIP_READERROR;
         }

         /* finding the variable associated with the column */
         getScenarioEntityName(name, getScenarioEntryCol(scenario, i), getScenarioStageNum(scenarioscip, scenario),
            getScenarioNum(scenarioscip, scenario));
         var = SCIPfindVar(scenarioscip, name);

         if( var == NULL )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", getScenarioEntryCol(scenario, i));
            var = SCIPfindVar(scenarioscip, name);
         }

         /* changing the coefficient for the variable */
         if( var == NULL )
         {
            SCIPerrorMessage("There is no variable <%s> in the current scenario.\n", name);
            return SCIP_READERROR;
         }
         else
         {
            SCIP_CALL( SCIPchgCoefLinear(scenarioscip, cons, var, getScenarioEntryValue(scenario, i)) );
         }
      }
   }

   return SCIP_OKAY;
}

/** removes the core variables and constriants for stage 2 and lower */
static
SCIP_RETCODE removeCoreVariablesAndConstraints(
   SCIP*                 scip                /**< the SCIP data structure */
   )
{
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   int nconss;
   int nvars;
   int numstages;
   int i;
   int j;
   SCIP_Bool deleted;

   assert(scip != NULL);

   numstages = SCIPtimGetNStages(scip);

   /* looping through all stages to remove the variables and constraints. The first stage is not removed as these are
    * part of the complete problem */
   for( i = 1; i < numstages; i++ )
   {
      conss = SCIPtimGetStageConss(scip, i);
      vars = SCIPtimGetStageVars(scip, i);
      nconss = SCIPtimGetStageNConss(scip, i);
      nvars = SCIPtimGetStageNVars(scip, i);

      /* removing constriants */
      for( j = 0; j < nconss; j++ )
      {
         if( !SCIPconsIsDeleted(conss[j]) )
            SCIP_CALL( SCIPdelCons(scip, conss[j]) );
      }

      /* removing variables */
      for( j = 0; j < nvars; j++ )
      {
         if( !SCIPvarIsDeleted(vars[j]) )
         {
            SCIP_CALL( SCIPdelVar(scip, vars[j], &deleted) );
            assert(deleted);
         }
      }
   }

   return SCIP_OKAY;
}


/* build the stochastic program completely as a MIP, i.e. no decomposition */
static
SCIP_RETCODE buildFullProblem(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   int i;

   assert(scip != NULL);
   assert(readerdata != NULL);

   /* adding all variables and constraints for stages below the first stage.
    * The first stage is covered by the original problem. */
   for( i = 0; i < getScenarioNChildren(readerdata->scenariotree); i++ )
      SCIP_CALL( addScenarioVarsAndConsToProb(scip, getScenarioChild(readerdata->scenariotree, i), FALSE) );

   /* removing the variable and constraints that were included as part of the core file */
   SCIP_CALL( removeCoreVariablesAndConstraints(scip) );

   return SCIP_OKAY;
}


/** builds the stochastic program using Benders' decomposition */
static
SCIP_RETCODE buildDecompProblem(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   SCIP_BENDERS* benders;
   int i;

   assert(scip != NULL);
   assert(readerdata != NULL);

   SCIP_CALL( createScenarioSubproblemArray(scip, readerdata->scenariotree) );

   /* activating the Benders' constraint handler. The two-phase method is activated by default. If the user desires not
    * to use the two-phase method, then the setting in cons_benderslp must be explicitly changed.
    */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/benders/active", TRUE) );

   setScenarioScip(readerdata->scenariotree, scip);

   /* adding all variables and constraints for stages below the first stage.
    * The first stage is covered by the original problem. */
   for( i = 0; i < getScenarioNChildren(readerdata->scenariotree); i++ )
   {
      SCIP_CALL( addScenarioVarsAndConsToProb(scip, getScenarioChild(readerdata->scenariotree, i), TRUE) );
      addScenarioSubproblem(readerdata->scenariotree, getScenarioScip(getScenarioChild(readerdata->scenariotree, i)));
   }

   /* creating the Benders' decomposition */
   SCIP_CALL( SCIPcreateBendersDefault(scip, getScenarioSubproblemArray(readerdata->scenariotree),
         getScenarioNChildren(readerdata->scenariotree)) );

   /* getting the default Benders' decomposition */
   benders = SCIPfindBenders(scip, "default");

   /* updating the lower bounds for the subproblems */
   for( i = 0; i < getScenarioNChildren(readerdata->scenariotree); i++ )
   {
      SCIPbendersUpdateSubproblemLowerbound(benders, i,
         getScenarioLowerbound(scip, getScenarioChild(readerdata->scenariotree, i)));
   }

   /* removing the variable and constraints that were included as part of the core file */
   SCIP_CALL( removeCoreVariablesAndConstraints(scip) );

   /* changing settings that are required for Benders' decomposition */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/trysol/freq", 1) );

   /* disabling aggregation since it can affect the mapping between the master and subproblem variables */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/donotaggr", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/donotmultaggr", TRUE) );

   return SCIP_OKAY;
}

/** Read the stochastic information of an SMPS file instance in "STO File Format". */
static
SCIP_RETCODE readSto(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of the input file */
   SCIP_READERDATA*      readerdata          /**< the reader data */
   )
{
   SCIP_FILE* fp;
   STOINPUT* stoi;
   SCIP_RETCODE retcode;
   SCIP_Bool error = TRUE;
   SCIP_Bool unsupported = FALSE;

   assert(scip != NULL);
   assert(filename != NULL);

   fp = SCIPfopen(filename, "r");
   if( fp == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   SCIP_CALL_FINALLY( stoinputCreate(scip, &stoi, fp), SCIPfclose(fp) );
   SCIP_CALL_TERMINATE( retcode, createReaderdata(scip, readerdata), TERMINATE );

   SCIP_CALL_TERMINATE( retcode, readStoch(scip, stoi), TERMINATE );

   /* checking for supported stochastic information types */
   if( stoinputStochInfoType(stoi) != STO_STOCHINFO_DISCRETE )
   {
         SCIPinfoMessage(scip, NULL, "\nSorry, currently only STO files with the stochastic information as DISCRETE are supported.\n\n");
         SCIPinfoMessage(scip, NULL, "NOTE: The problem provided by the COR file is loaded without stochastic information.\n\n");
         unsupported = TRUE;
   }
   else
   {
      if( stoinputSection(stoi) == STO_BLOCKS )
      {
         /* coverity[tainted_data] */
         SCIP_CALL_TERMINATE( retcode, readBlocks(stoi, scip, readerdata), TERMINATE );
      }

      if( stoinputSection(stoi) == STO_SCENARIOS )
      {
         /* if there are more than two stages, then the sto file is not read. */
         if( SCIPtimGetNStages(scip) > 2 )
         {
            SCIPinfoMessage(scip, NULL, "\nThe scenarios for the stochastic programs are defined in <%s> as SCENARIOS\n", filename);
            SCIPinfoMessage(scip, NULL, "Sorry, currently only two-stage stochastic programs are supported when scenarios are defined as SCENARIOS.\n\n");
            SCIPinfoMessage(scip, NULL, "NOTE: The problem provided by the COR file is loaded without stochastic information.\n\n");
            unsupported = TRUE;
         }
         else
         {
            SCIP_CALL_TERMINATE( retcode, readScenarios(stoi, scip, readerdata), TERMINATE );
         }
      }

      if( stoinputSection(stoi) == STO_INDEP )
      {
         SCIP_CALL_TERMINATE( retcode, readIndep(stoi, scip, readerdata), TERMINATE );
      }
   }

   if( !unsupported && stoinputSection(stoi) != STO_ENDATA )
      stoinputSyntaxerror(stoi);

   error = stoinputHasError(stoi);

   if( !error && !unsupported )
   {
      if( readerdata->usebenders )
      {
         SCIP_CALL_TERMINATE( retcode, buildDecompProblem(scip, readerdata), TERMINATE );
      }
      else
      {
         SCIP_CALL_TERMINATE( retcode, buildFullProblem(scip, readerdata), TERMINATE );
      }
   }

/* cppcheck-suppress unusedLabel */
TERMINATE:
   stoinputFree(scip, &stoi);
   SCIPfclose(fp);

   if( error || retcode != SCIP_OKAY )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySto)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderSto(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeSto)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( freeReaderdata(scip, readerdata) );

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSto)
{  /*lint --e{715}*/
   SCIP_READER* correader;
   SCIP_READER* timreader;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   correader = SCIPfindReader(scip, "correader");
   timreader = SCIPfindReader(scip, "timreader");

   if( correader == NULL )
   {
      SCIPwarningMessage(scip, "It is necessary to include the \"cor\" reader\n");
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( timreader == NULL )
   {
      SCIPwarningMessage(scip, "It is necessary to include the \"tim\" reader\n");
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* checking whether the cor file has been read */
   if( !SCIPcorHasRead(correader) )
   {
      SCIPwarningMessage(scip, "The core file must be read before the time and stochastic files.\n");
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* checking whether the tim file has been read */
   if( !SCIPtimHasRead(timreader) )
   {
      SCIPwarningMessage(scip, "The time file must be read before the stochastic files.\n");
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPreadSto(scip, filename, result) );

   return SCIP_OKAY;
}

/*
 * sto file reader specific interface methods
 */

/** includes the sto file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSto(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );
   readerdata->scenariotree = NULL;
   readerdata->numscenarios = 0;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySto) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeSto) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSto) );

   /* add decomposition parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/usebenders",
         "should Benders' decomposition be used?",
         &readerdata->usebenders, FALSE, DEFAULT_USEBENDERS, NULL, NULL) );

   return SCIP_OKAY;
}


/** reads the stochastic information for a stochastic program that is in SMPS format */
SCIP_RETCODE SCIPreadSto(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(result != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);
   readerdata = SCIPreaderGetData(reader);

   retcode = readSto(scip, filename, readerdata);

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   if( retcode == SCIP_NOFILE || retcode == SCIP_READERROR )
      return retcode;

   SCIP_CALL( retcode );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** returns the total number of scenarios added to the problem */
int SCIPstoGetNScenarios(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->numscenarios;
}
