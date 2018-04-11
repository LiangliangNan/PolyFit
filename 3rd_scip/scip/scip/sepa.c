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

/**@file   sepa.c
 * @brief  methods and datastructures for separators
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/sepastore.h"
#include "scip/scip.h"
#include "scip/sepa.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_sepa.h"


/** compares two separators w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPsepaComp)
{  /*lint --e{715}*/
   return ((SCIP_SEPA*)elem2)->priority - ((SCIP_SEPA*)elem1)->priority;
}

/** comparison method for sorting separators w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPsepaCompName)
{
   return strcmp(SCIPsepaGetName((SCIP_SEPA*)elem1), SCIPsepaGetName((SCIP_SEPA*)elem2));
}

/** method to call, when the priority of a separator was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdSepaPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSepaPriority() to mark the sepas unsorted */
   SCIP_CALL( SCIPsetSepaPriority(scip, (SCIP_SEPA*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given separator to a new scip */
SCIP_RETCODE SCIPsepaCopyInclude(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( sepa->sepacopy != NULL )
   {
      SCIPsetDebugMsg(set, "including separator %s in subscip %p\n", SCIPsepaGetName(sepa), (void*)set->scip);
      SCIP_CALL( sepa->sepacopy(set->scip, sepa) );
   }
   return SCIP_OKAY;
}

/** creates a separator */
SCIP_RETCODE SCIPsepaCreate(
   SCIP_SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPACOPY    ((*sepacopy)),      /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   SCIP_DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol)),   /**< solving process initialization method of separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol)),   /**< solving process deinitialization method of separator */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp)),    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol)),   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(sepa != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(0.0 <= maxbounddist && maxbounddist <= 1.0);
   assert(sepaexeclp != NULL || sepaexecsol != NULL);

   SCIP_ALLOC( BMSallocMemory(sepa) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*sepa)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*sepa)->desc, desc, strlen(desc)+1) );
   (*sepa)->priority = priority;
   (*sepa)->freq = freq;
   (*sepa)->maxbounddist = maxbounddist;
   (*sepa)->usessubscip = usessubscip;
   (*sepa)->sepacopy = sepacopy;
   (*sepa)->sepafree = sepafree;
   (*sepa)->sepainit = sepainit;
   (*sepa)->sepaexit = sepaexit;
   (*sepa)->sepainitsol = sepainitsol;
   (*sepa)->sepaexitsol = sepaexitsol;
   (*sepa)->sepaexeclp = sepaexeclp;
   (*sepa)->sepaexecsol = sepaexecsol;
   (*sepa)->sepadata = sepadata;
   SCIP_CALL( SCIPclockCreate(&(*sepa)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*sepa)->sepaclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*sepa)->lastsepanode = -1;
   (*sepa)->ncalls = 0;
   (*sepa)->ncutoffs = 0;
   (*sepa)->ncutsfound = 0;
   (*sepa)->ncutsapplied = 0;
   (*sepa)->nconssfound = 0;
   (*sepa)->ndomredsfound = 0;
   (*sepa)->ncallsatnode = 0;
   (*sepa)->ncutsfoundatnode = 0;
   (*sepa)->lpwasdelayed = FALSE;
   (*sepa)->solwasdelayed = FALSE;
   (*sepa)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of separator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*sepa)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4,
         paramChgdSepaPriority, (SCIP_PARAMDATA*)(*sepa)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/freq", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "frequency for calling separator <%s> (-1: never, 0: only in root node)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*sepa)->freq, FALSE, freq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/maxbounddist", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying separator <%s> (0.0: only on current best node, 1.0: on all nodes)",
      name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*sepa)->maxbounddist, TRUE, maxbounddist, 0.0, 1.0, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/delay", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should separator be delayed, if other separators found cuts?",
         &(*sepa)->delay, TRUE, delay, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "separating/%s/expbackoff", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "base for exponential increase of frequency at which separator <%s> is called (1: call at each multiple of frequency)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*sepa)->expbackoff, TRUE, 4, 1, 100, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of separator */
SCIP_RETCODE SCIPsepaFree(
   SCIP_SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(*sepa != NULL);
   assert(!(*sepa)->initialized);
   assert(set != NULL);

   /* call destructor of separator */
   if( (*sepa)->sepafree != NULL )
   {
      SCIP_CALL( (*sepa)->sepafree(set->scip, *sepa) );
   }

   SCIPclockFree(&(*sepa)->sepaclock);
   SCIPclockFree(&(*sepa)->setuptime);
   BMSfreeMemoryArray(&(*sepa)->name);
   BMSfreeMemoryArray(&(*sepa)->desc);
   BMSfreeMemory(sepa);

   return SCIP_OKAY;
}

/** initializes separator */
SCIP_RETCODE SCIPsepaInit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   if( sepa->initialized )
   {
      SCIPerrorMessage("separator <%s> already initialized\n", sepa->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(sepa->setuptime);
      SCIPclockReset(sepa->sepaclock);

      sepa->lastsepanode = -1;
      sepa->ncalls = 0;
      sepa->ncutoffs = 0;
      sepa->ncutsfound = 0;
      sepa->ncutsapplied = 0;
      sepa->nconssfound = 0;
      sepa->ndomredsfound = 0;
      sepa->ncallsatnode = 0;
      sepa->ncutsfoundatnode = 0;
      sepa->lpwasdelayed = FALSE;
      sepa->solwasdelayed = FALSE;
   }

   if( sepa->sepainit != NULL )
   {
      /* start timing */
      SCIPclockStart(sepa->setuptime, set);

      SCIP_CALL( sepa->sepainit(set->scip, sepa) );

      /* stop timing */
      SCIPclockStop(sepa->setuptime, set);
   }
   sepa->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of separator */
SCIP_RETCODE SCIPsepaExit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   if( !sepa->initialized )
   {
      SCIPerrorMessage("separator <%s> not initialized\n", sepa->name);
      return SCIP_INVALIDCALL;
   }

   if( sepa->sepaexit != NULL )
   {
      /* start timing */
      SCIPclockStart(sepa->setuptime, set);

      SCIP_CALL( sepa->sepaexit(set->scip, sepa) );

      /* stop timing */
      SCIPclockStop(sepa->setuptime, set);
   }
   sepa->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs separator that the branch and bound process is being started */
SCIP_RETCODE SCIPsepaInitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   sepa->lpwasdelayed = FALSE;
   sepa->solwasdelayed = FALSE;

   /* call solving process initialization method of separator */
   if( sepa->sepainitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(sepa->setuptime, set);

      SCIP_CALL( sepa->sepainitsol(set->scip, sepa) );

      /* stop timing */
      SCIPclockStop(sepa->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs separator that the branch and bound process data is being freed */
SCIP_RETCODE SCIPsepaExitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of separator */
   if( sepa->sepaexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(sepa->setuptime, set);

      SCIP_CALL( sepa->sepaexitsol(set->scip, sepa) );

      /* stop timing */
      SCIPclockStop(sepa->setuptime, set);
   }

   return SCIP_OKAY;
}

/** calls LP separation method of separator */
SCIP_RETCODE SCIPsepaExecLP(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   depth,              /**< depth of current node */
   SCIP_Real             bounddist,          /**< current relative distance of local dual bound to global dual bound */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(sepa != NULL);
   assert(sepa->freq >= -1);
   assert(0.0 <= sepa->maxbounddist && sepa->maxbounddist <= 1.0);
   assert(0.0 <= bounddist && bounddist <= 1.0);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( sepa->sepaexeclp != NULL && SCIPsetIsLE(set, bounddist, sepa->maxbounddist) &&
       ( (depth == 0 && sepa->freq != -1) ||
         (sepa->freq > 0 && depth % sepa->freq == 0 &&
            (sepa->expbackoff == 1 || SCIPsetIsIntegral(set, LOG2(depth * (1.0 / sepa->freq)) / LOG2((SCIP_Real)sepa->expbackoff)))) ||
         sepa->lpwasdelayed )
     )
   {
      if( (!sepa->delay && !sepa->lpwasdelayed) || execdelayed )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         int oldncuts;
         int oldnactiveconss;
         int ncutsfound;

         SCIPsetDebugMsg(set, "executing separator <%s> on LP solution\n", sepa->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         oldncuts = SCIPsepastoreGetNCuts(sepastore);
         oldnactiveconss = stat->nactiveconss;

         /* reset the statistics for current node */
         if( sepa->lastsepanode != stat->ntotalnodes )
         {
            sepa->ncallsatnode = 0;
            sepa->ncutsfoundatnode = 0;
         }

         /* start timing */
         SCIPclockStart(sepa->sepaclock, set);

         /* call external separation method */
         SCIP_CALL( sepa->sepaexeclp(set->scip, sepa, result, allowlocal) );

         /* stop timing */
         SCIPclockStop(sepa->sepaclock, set);

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
         {
            sepa->ncalls++;
            sepa->ncallsatnode++;
            sepa->lastsepanode = stat->ntotalnodes;
         }
         if( *result == SCIP_CUTOFF )
            sepa->ncutoffs++;
         ncutsfound = SCIPsepastoreGetNCuts(sepastore) - oldncuts;
         sepa->ncutsfound += ncutsfound;
         sepa->ncutsfoundatnode += ncutsfound;
         sepa->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         sepa->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         sepa->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_NEWROUND
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            SCIPerrorMessage("execution method of separator <%s> returned invalid result <%d>\n",
               sepa->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "separator <%s> was delayed\n", sepa->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether separator was delayed */
      sepa->lpwasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** calls primal solution separation method of separator */
SCIP_RETCODE SCIPsepaExecSol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             allowlocal,         /**< should the separator allow local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(sepa != NULL);
   assert(sepa->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( sepa->sepaexecsol != NULL &&
       ( (depth == 0 && sepa->freq != -1) ||
         (sepa->freq > 0 && depth % sepa->freq == 0 &&
            (sepa->expbackoff == 1 || SCIPsetIsIntegral(set, LOG2(depth * (1.0 / sepa->freq) / LOG2((SCIP_Real)sepa->expbackoff))))) ||
         sepa->solwasdelayed )
     )
   {
      if( (!sepa->delay && !sepa->solwasdelayed) || execdelayed )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         int oldncuts;
         int oldnactiveconss;
         int ncutsfound;

         SCIPsetDebugMsg(set, "executing separator <%s> on solution %p\n", sepa->name, (void*)sol);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         oldncuts = SCIPsepastoreGetNCuts(sepastore);
         oldnactiveconss = stat->nactiveconss;

         /* reset the statistics for current node */
         if( sepa->lastsepanode != stat->ntotalnodes )
         {
            sepa->ncallsatnode = 0;
            sepa->ncutsfoundatnode = 0;
         }

         /* start timing */
         SCIPclockStart(sepa->sepaclock, set);

         /* call external separation method */
         SCIP_CALL( sepa->sepaexecsol(set->scip, sepa, sol, result, allowlocal) );

         /* stop timing */
         SCIPclockStop(sepa->sepaclock, set);

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
         {
            sepa->ncalls++;
            sepa->ncallsatnode++;
            sepa->lastsepanode = stat->ntotalnodes;
         }
         if( *result == SCIP_CUTOFF )
            sepa->ncutoffs++;
         ncutsfound = SCIPsepastoreGetNCuts(sepastore) - oldncuts;
         sepa->ncutsfound += ncutsfound;
         sepa->ncutsfoundatnode += ncutsfound;
         sepa->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         sepa->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         sepa->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_NEWROUND
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            SCIPerrorMessage("execution method of separator <%s> returned invalid result <%d>\n",
               sepa->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "separator <%s> was delayed\n", sepa->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether separator was delayed */
      sepa->solwasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** gets user data of separator */
SCIP_SEPADATA* SCIPsepaGetData(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->sepadata;
}

/** sets user data of separator; user has to free old data in advance! */
void SCIPsepaSetData(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata            /**< new separator user data */
   )
{
   assert(sepa != NULL);

   sepa->sepadata = sepadata;
}

/* new callback/method setter methods */

/** sets copy method of separator */
void SCIPsepaSetCopy(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPACOPY    ((*sepacopy))       /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(sepa != NULL);

   sepa->sepacopy = sepacopy;
}

/** sets destructor method of separator */
void SCIPsepaSetFree(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAFREE    ((*sepafree))       /**< destructor of separator */
   )
{
   assert(sepa != NULL);

   sepa->sepafree = sepafree;
}

/** sets initialization method of separator */
void SCIPsepaSetInit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINIT    ((*sepainit))       /**< initialize separator */
   )
{
   assert(sepa != NULL);

   sepa->sepainit = sepainit;
}

/** sets deinitialization method of separator */
void SCIPsepaSetExit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit))       /**< deinitialize separator */
   )
{
   assert(sepa != NULL);

   sepa->sepaexit = sepaexit;
}

/** sets solving process initialization method of separator */
void SCIPsepaSetInitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol))    /**< solving process initialization method of separator */
   )
{
   assert(sepa != NULL);

   sepa->sepainitsol = sepainitsol;
}

/** sets solving process deinitialization method of separator */
void SCIPsepaSetExitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol))    /**< solving process deinitialization method of separator */
   )
{
   assert(sepa != NULL);

   sepa->sepaexitsol = sepaexitsol;
}

/** gets name of separator */
const char* SCIPsepaGetName(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->name;
}

/** gets description of separator */
const char* SCIPsepaGetDesc(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->desc;
}

/** gets priority of separator */
int SCIPsepaGetPriority(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->priority;
}

/** sets priority of separator */
void SCIPsepaSetPriority(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the separator */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   sepa->priority = priority;
   set->sepassorted = FALSE;
}

/** gets frequency of separator */
int SCIPsepaGetFreq(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->freq;
}

/** sets frequency of separator */
void SCIPsepaSetFreq(
   SCIP_SEPA*            sepa,               /**< separator */
   int                   freq                /**< new frequency of separator */
   )
{
   assert(sepa != NULL);

   sepa->freq = freq;
}

/** get maximal bound distance at which the separator is called */
SCIP_Real SCIPsepaGetMaxbounddist(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->maxbounddist;
}

/** does the separator use a secondary SCIP instance? */
SCIP_Bool SCIPsepaUsesSubscip(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->usessubscip;
}

/** enables or disables all clocks of \p sepa, depending on the value of the flag */
void SCIPsepaEnableOrDisableClocks(
   SCIP_SEPA*            sepa,               /**< the separator for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the separator be enabled? */
   )
{
   assert(sepa != NULL);

   SCIPclockEnableOrDisable(sepa->setuptime, enable);
   SCIPclockEnableOrDisable(sepa->sepaclock, enable);
}

/** gets time in seconds used in this separator for setting up for next stages */
SCIP_Real SCIPsepaGetSetupTime(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return SCIPclockGetTime(sepa->setuptime);
}

/** gets time in seconds used in this separator */
SCIP_Real SCIPsepaGetTime(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return SCIPclockGetTime(sepa->sepaclock);
}

/** gets the total number of times, the separator was called */
SCIP_Longint SCIPsepaGetNCalls(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncalls;
}

/** gets the number of times, the separator was called at the current node */
int SCIPsepaGetNCallsAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncallsatnode;
}

/** gets total number of times, the separator detected a cutoff */
SCIP_Longint SCIPsepaGetNCutoffs(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutoffs;
}

/** gets the total number of cutting planes found by this separator */
SCIP_Longint SCIPsepaGetNCutsFound(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsfound;
}

/** gets the total number of cutting planes applied to lp */
SCIP_Longint SCIPsepaGetNCutsApplied(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsapplied;
}

/** increase count of applied cuts */
void SCIPsepaIncNAppliedCuts(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert( sepa != NULL );

   ++sepa->ncutsapplied;
}

/** increase count of found cuts */
void SCIPsepaIncNCutsFound(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert( sepa != NULL );

   ++sepa->ncutsfound;
}

/** increase count of found cuts at current node */
void SCIPsepaIncNCutsFoundAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert( sepa != NULL );

   ++sepa->ncutsfoundatnode;
}

/** gets the number of cutting planes found by this separator at the current node */
SCIP_Longint SCIPsepaGetNCutsFoundAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsfoundatnode;
}

/** gets total number of additional constraints added by this separator */
SCIP_Longint SCIPsepaGetNConssFound(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->nconssfound;
}

/** gets total number of domain reductions found by this separator */
SCIP_Longint SCIPsepaGetNDomredsFound(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ndomredsfound;
}

/** should separator be delayed, if other separators found cuts? */
SCIP_Bool SCIPsepaIsDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->delay;
}

/** was separation of the LP solution delayed at the last call? */
SCIP_Bool SCIPsepaWasLPDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->lpwasdelayed;
}

/** was separation of the primal solution delayed at the last call? */
SCIP_Bool SCIPsepaWasSolDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->solwasdelayed;
}

/** is separator initialized? */
SCIP_Bool SCIPsepaIsInitialized(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->initialized;
}
