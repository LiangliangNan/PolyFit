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

/**@file   presol_gateextraction.c
 * @brief  gateextraction presolver
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_gateextraction.h"
#include "scip/cons_setppc.h"
#include "scip/cons_logicor.h"
#include "scip/cons_and.h"


#define PRESOL_NAME            "gateextraction"
#define PRESOL_DESC            "presolver extracting gate(and)-constraints"
#define PRESOL_PRIORITY         1000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define HASHSIZE_LOGICORCONS     500 /**< minimal size of hash table in logicor constraint tables */
#define HASHSIZE_SETPPCCONS      500 /**< minimal size of hash table in setppc constraint tables */

#define DEFAULT_ONLYSETPART       FALSE  /**< should only set-partitioning constraints be extracted and no and-constraints */
#define DEFAULT_SEARCHEQUATIONS    TRUE  /**< should we try to extract set-partitioning constraint out of one logicor
					  *   and one corresponding set-packing constraint
					  */
#define DEFAULT_SORTING               1  /**< order logicor contraints to extract big-gates before smaller ones (-1), do
					  *   not order them (0) or order them to extract smaller gates at first (1)
					  */


/* This presolver tries to extract gate-constraints meaning and-constraints and set-partitioning constraints (and could
 * be expanded to find xor-constraints too). This is done by detecting linearizations or systems of inequalities which
 * form an and-constraint or a set-partitioning constraint. An example:
 *
 * we have a logicor constraint of the form:                x + y + z >= 1
 *
 * and we also have the following set-packing constraints: (x + y <= 1 and x + z <= 1) <=> (~x + ~y >= 1 and ~x + ~z >= 1)
 *
 * - these three constraints form an and-constraint:        x = ~y * ~z (x = AND(~y,~z))
 *
 * if an additional set-packing constraint exists:          y + z <= 1
 *
 * - these four constraints form a set-partitioning cons.:  x + y + z = 1
 *
 * some information can be found:
 *
 *  http://www.cs.ubc.ca/~hutter/earg/papers07/cnf-structure.pdf
 *  http://www.cadence.com/cn/cadence/cadence_labs/Documents/niklas_SAT_2005_Effective.pdf
 *
 * We also do some check for logicor and set-packing/-partitioning constraint with the same variables to upgrade these
 * both constraints into one. For example:
 *
 *  x + y + z >= 1 and x + y + z <= 1 form x + y + z = 1
 *
 */


/*
 * Data structures
 */


/** data object to compare constraint easier */
struct HashData
{
   SCIP_CONS*            cons;               /**< pointer the the corresponding constraint */
   SCIP_VAR**            vars;               /**< constraint variables used for hash comparison */
   int                   nvars;              /**< number of variables */
};
typedef struct HashData HASHDATA;


/** presolver data */
struct SCIP_PresolData
{
   HASHDATA*             setppchashdatas;    /**< setppc-hashdata storage */
   SCIP_HASHTABLE*       hashdatatable;      /**< setppc-hashdata hashtable for usable setppc constraints */
   SCIP_HASHTABLE*       setppchashtable;    /**< setppc hashtable for usable setppc constraints */
   SCIP_HASHTABLE*       logicorhashtable;   /**< logicor hashtable for usable logicor constraints */
   SCIP_CONS**           usefullogicor;      /**< array for usable logicors */
   int                   nusefullogicor;     /**< number of usable logicors */
   int                   susefullogicor;     /**< size of array for usable logicor constraints */
   int                   nsetppchashdatas;   /**< number of setppchashdata elements added to the hashtable */
   int                   ssetppchashdatas;   /**< size of setppchashdata elements added to the hashtable */
   int                   ngates;             /**< number of found gates in presolving */
   int                   firstchangedlogicor;/**< position of the first new/changed logicor constraint in the
                                              *   usefullogicor array
                                              */
   int                   maxnvarslogicor;    /**< maximal number of variables a logicor constraint has */
   int                   sorting;            /**< integer parameter how to sort logicor constraints for extracting gates */
   SCIP_Bool             usefulsetppcexist;  /**< did we find usable set-packing constraints for gate extraction */
   SCIP_Bool             usefullogicorexist; /**< did we find usable logicor constraints for gate extraction */
   SCIP_Bool             newsetppchashdatas; /**< flag indicating whether we found new set-packing constraint with two
                                              *   variables since the last presolving round
                                              */
   SCIP_Bool             initialized;        /**< was data alredy be initialized */
   SCIP_Bool             onlysetpart;        /**< boolean parameter whetehr we only want to extract linear gates */
   SCIP_Bool             searchequations;    /**< boolean parameter whetehr we want to search for equations arising from
                                              *   logicor and setppc constraints
                                              */
};


/*
 * Local methods
 */


/** returns TRUE iff both keys are equal; two constraints are equal if they have the same pointer */
static
SCIP_DECL_HASHKEYEQ(hashdataKeyEqCons)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   HASHDATA* hashdata1;
   HASHDATA* hashdata2;
   int v;

   hashdata1 = (HASHDATA*)key1;
   hashdata2 = (HASHDATA*)key2;
#ifndef NDEBUG
   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   /* check data structure */
   assert(hashdata1->nvars == 2);
   assert(hashdata2->nvars == 2);
   /* at least one data object needs to be have a real set packing constraint */
   /* TODO why does this assert fail on one instance when problem is freed
    * using the new hashing: assert(hashdata1->cons != NULL || hashdata2->cons != NULL);
    */

   for( v = 1; v >= 0; --v )
   {
      /* tests if variables are equal */
      if( hashdata1->vars[v] != hashdata2->vars[v] )
         return FALSE;

      assert(SCIPvarCompare(hashdata1->vars[v], hashdata2->vars[v]) == 0);
   }

   /* a hashdata object is only equal if it has the same constraint pointer, or one has no constraint pointer, latter
    * means that this hashdata object is derived from a logicor constraint
    */
   if( hashdata1->cons == NULL || hashdata2->cons == NULL || hashdata1->cons == hashdata2->cons )
      return TRUE;
   else
      return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashdataKeyValCons)
{  /*lint --e{715}*/
   HASHDATA* hashdata;
   unsigned int hashval;

   hashdata = (HASHDATA*)key;
   assert(hashdata != NULL);
   assert(hashdata->vars != NULL);
   assert(hashdata->nvars == 2);

   /* if we have only two variables we store at each 16 bits of the hash value the index of a variable */
   hashval = ((unsigned int)SCIPvarGetIndex(hashdata->vars[1]) << 16) + SCIPvarGetIndex(hashdata->vars[0]); /*lint !e701*/

   return hashval;
}


/** returns TRUE iff both keys are equal; two constraints are equal if they have the same pointer */
static
SCIP_DECL_HASHKEYEQ(setppcHashdataKeyEqCons)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   HASHDATA* hashdata1;
   HASHDATA* hashdata2;
   int v;

   hashdata1 = (HASHDATA*)key1;
   hashdata2 = (HASHDATA*)key2;
#ifndef NDEBUG
   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   /* check data structure */
   assert(hashdata1->nvars >= 2);
   assert(hashdata2->nvars >= 2);
   /* at least one data object needs to be have a real set-packing/partitioning constraint */
   assert(hashdata1->cons != NULL || hashdata2->cons != NULL);

   if( hashdata1->nvars != hashdata2->nvars )
      return FALSE;

   for( v = hashdata1->nvars - 1; v >= 0; --v )
   {
      /* tests if variables are equal */
      if( hashdata1->vars[v] != hashdata2->vars[v] )
         return FALSE;

      assert(SCIPvarCompare(hashdata1->vars[v], hashdata2->vars[v]) == 0);
   }

   /* a hashdata object is only equal if it has the same constraint pointer, or one has no constraint pointer, latter
    * means that this hashdata object is derived from a logicor constraint
    */
   if( hashdata1->cons == NULL || hashdata2->cons == NULL || hashdata1->cons == hashdata2->cons )
      return TRUE;
   else
      return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(setppcHashdataKeyValCons)
{  /*lint --e{715}*/
   HASHDATA* hashdata;

   hashdata = (HASHDATA*)key;
   assert(hashdata != NULL);
   assert(hashdata->vars != NULL);
   assert(hashdata->nvars >= 2);

   return SCIPhashTwo(SCIPcombineTwoInt(hashdata->nvars, SCIPvarGetIndex(hashdata->vars[0])), \
                      SCIPcombineTwoInt(SCIPvarGetIndex(hashdata->vars[hashdata->nvars/2]), \
                                        SCIPvarGetIndex(hashdata->vars[hashdata->nvars-1])));
}

/** initialize gateextraction presolver data */
static
void presoldataInit(
   SCIP_PRESOLDATA*      presoldata          /**< data object of presolver */
   )
{
   assert(presoldata != NULL);

   presoldata->usefullogicor = NULL;
   presoldata->nusefullogicor = 0;
   presoldata->susefullogicor = 0;
   presoldata->firstchangedlogicor = -1;
   presoldata->maxnvarslogicor = 0;;
   presoldata->nsetppchashdatas = 0;
   presoldata->ssetppchashdatas = 0;
   presoldata->ngates = 0;
   presoldata->usefulsetppcexist = FALSE;
   presoldata->usefullogicorexist = FALSE;
   presoldata->newsetppchashdatas = FALSE;
   presoldata->initialized = FALSE;

   presoldata->hashdatatable = NULL;
   presoldata->setppchashtable = NULL;
   presoldata->logicorhashtable = NULL;
}

/** initialize gateextraction hashtables */
static
SCIP_RETCODE presoldataInitHashtables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< data object of presolver */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);

   assert(presoldata->nusefullogicor == 0);
   assert(presoldata->susefullogicor == 0);
   assert(presoldata->nsetppchashdatas == 0);
   assert(presoldata->ssetppchashdatas == 0);
   assert(presoldata->firstchangedlogicor == -1);
   assert(presoldata->ngates == 0);
   assert(presoldata->usefullogicorexist == FALSE);
   assert(presoldata->usefulsetppcexist == FALSE);
   assert(presoldata->newsetppchashdatas == FALSE);
   assert(presoldata->initialized == FALSE);

   assert(presoldata->hashdatatable == NULL);
   assert(presoldata->setppchashtable == NULL);
   assert(presoldata->logicorhashtable == NULL);

   /* create hashtables */
   SCIP_CALL( SCIPhashtableCreate(&(presoldata->hashdatatable), SCIPblkmem(scip), HASHSIZE_SETPPCCONS,
                                  SCIPhashGetKeyStandard, hashdataKeyEqCons, hashdataKeyValCons, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&(presoldata->setppchashtable), SCIPblkmem(scip), HASHSIZE_SETPPCCONS,
                                  SCIPhashGetKeyStandard, SCIPhashKeyEqPtr, SCIPhashKeyValPtr, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&(presoldata->logicorhashtable), SCIPblkmem(scip), HASHSIZE_LOGICORCONS,
                                  SCIPhashGetKeyStandard, SCIPhashKeyEqPtr, SCIPhashKeyValPtr, (void*) scip) );

   return SCIP_OKAY;
}


/** create useful set-packing information by adding new set-packing constraints with two variables */
static
SCIP_RETCODE createPresoldata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< data object of presolver */
   SCIP_CONS**           setppcs,            /**< active setppc constraints */
   int                   nsetppcs,           /**< number of active setppc constraints */
   SCIP_CONS**           logicors,           /**< active logicor constraints */
   int                   nlogicors           /**< number of active logicor constraints */
   )
{
   SCIP_CONS** usefulconss;
   int nusefulconss = 0;
   int size;
   int c;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(setppcs != NULL);
   assert(nsetppcs > 0);
   assert(logicors != NULL);
   assert(nlogicors > 0);
   assert(presoldata->setppchashtable != NULL);
   assert(presoldata->logicorhashtable != NULL);

   presoldata->initialized = TRUE;

   size = MAX(nsetppcs, nlogicors);

   /* temporary memory for collecting set-packing constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &usefulconss, size) );

   if( !presoldata->usefulsetppcexist )
   {
      /* find set-packing constraints with exactly two variables */
      for( c = 0; c < nsetppcs; ++c )
      {
         assert(SCIPconsIsActive(setppcs[c]));

         if( SCIPgetTypeSetppc(scip, setppcs[c]) == SCIP_SETPPCTYPE_PACKING && SCIPgetNVarsSetppc(scip, setppcs[c]) == 2 && !SCIPconsIsModifiable(setppcs[c]) )
         {
            /* insert new element in hashtable */
            SCIP_CALL( SCIPhashtableInsert(presoldata->setppchashtable, (void*) setppcs[c]) );

            usefulconss[nusefulconss] = setppcs[c];
            ++nusefulconss;
         }
      }

      /* add usefulconss constraints to hashdata elements */
      if( nusefulconss > 0 )
      {
         SCIP_Bool negated[2];
         int h;

         presoldata->usefulsetppcexist = TRUE;
         presoldata->ssetppchashdatas = nusefulconss;

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->setppchashdatas), nusefulconss) );

         h = 0;
         for( c = 0; c < nusefulconss; ++c )
         {
            SCIP_VAR** setppcvars = SCIPgetVarsSetppc(scip, usefulconss[c]);
            assert(SCIPconsIsActive(usefulconss[c]));
            assert(SCIPgetNVarsSetppc(scip, usefulconss[c]) == 2);

            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(presoldata->setppchashdatas[h].vars), setppcvars, 2) );

            SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[h].vars[0], &(presoldata->setppchashdatas[h].vars[0]), &(negated[0])) );
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[h].vars[1], &(presoldata->setppchashdatas[h].vars[1]), &(negated[1])) );

            if( SCIPvarGetStatus(presoldata->setppchashdatas[h].vars[0]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[h].vars[0]) == SCIP_VARSTATUS_MULTAGGR
                  || SCIPvarGetStatus(presoldata->setppchashdatas[h].vars[1]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[h].vars[1]) == SCIP_VARSTATUS_MULTAGGR )
            {
               SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[h].vars), 2);
               continue;
            }

            presoldata->setppchashdatas[h].nvars = 2;

            /* capture variables */
            SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[h].vars[0]) );
            SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[h].vars[1]) );

            /* order the variables after their index */
            if( SCIPvarGetIndex(presoldata->setppchashdatas[h].vars[0]) > SCIPvarGetIndex(presoldata->setppchashdatas[h].vars[1]) )
            {
               SCIP_VAR* tmp = presoldata->setppchashdatas[h].vars[0];
               presoldata->setppchashdatas[h].vars[0] = presoldata->setppchashdatas[h].vars[1];
               presoldata->setppchashdatas[h].vars[1] = tmp;
            }

            presoldata->setppchashdatas[h].cons = usefulconss[c];

            SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[h]) );
            SCIP_CALL( SCIPcaptureCons(scip, usefulconss[c]) );

            ++h;
         }
         presoldata->nsetppchashdatas = h;

         if( presoldata->nsetppchashdatas > 0 )
            presoldata->newsetppchashdatas = TRUE;
      }
   }

   nusefulconss = 0;

   if( !presoldata->usefullogicorexist )
   {
      /* capture all logicor constraints */
      for( c = 0; c < nlogicors; ++c )
      {
         assert(SCIPconsIsActive(logicors[c]));

         if( !SCIPconsIsModifiable(logicors[c]) && SCIPgetNVarsLogicor(scip, logicors[c]) >= 3 )
         {
            /* insert new element in hashtable */
            SCIP_CALL( SCIPhashtableInsert(presoldata->logicorhashtable, (void*) logicors[c]) );
            SCIP_CALL( SCIPcaptureCons(scip, logicors[c]) );

            usefulconss[nusefulconss] = logicors[c];
            ++nusefulconss;

            /* update maximal entries in a logicor constraint */
            if( presoldata->maxnvarslogicor < SCIPgetNVarsLogicor(scip, logicors[c]) )
               presoldata->maxnvarslogicor = SCIPgetNVarsLogicor(scip, logicors[c]);
         }
      }

      /* no usefulconss constraints */
      if( nusefulconss > 0 )
      {
         presoldata->firstchangedlogicor = 0;
         presoldata->usefullogicorexist = TRUE;
         presoldata->susefullogicor = nusefulconss;
         presoldata->nusefullogicor = nusefulconss;
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &presoldata->usefullogicor, usefulconss, presoldata->susefullogicor) );
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &usefulconss);

   return SCIP_OKAY;
}


/** remove old setppchashdatas objects, so that the allocated memory will stay low */
static
SCIP_RETCODE cleanupHashDatas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< data object of presolver */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);

   if( presoldata->usefulsetppcexist )
   {
      int c;

      assert(presoldata->setppchashdatas != NULL || presoldata->nsetppchashdatas == 0);

      for( c = presoldata->nsetppchashdatas - 1; c >= 0; --c )
      {
         SCIP_Bool removeentry = FALSE;

         assert(presoldata->setppchashdatas[c].cons != NULL);

         if( SCIPconsIsDeleted(presoldata->setppchashdatas[c].cons) || SCIPconsIsModifiable(presoldata->setppchashdatas[c].cons)
               || SCIPgetTypeSetppc(scip, presoldata->setppchashdatas[c].cons) != SCIP_SETPPCTYPE_PACKING || SCIPgetNVarsSetppc(scip, presoldata->setppchashdatas[c].cons) != 2 )
         {
            removeentry = TRUE;
         }
         else
         {
            SCIP_VAR* vars[2];
            SCIP_Bool negated[2];

            SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[c].vars[0], &(vars[0]), &(negated[0])) );
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[c].vars[1], &(vars[1]), &(negated[1])) );

            if( SCIPvarGetStatus(vars[0]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(vars[0]) == SCIP_VARSTATUS_MULTAGGR
                  || SCIPvarGetStatus(vars[1]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(vars[1]) == SCIP_VARSTATUS_MULTAGGR
                  || presoldata->setppchashdatas[c].vars[0] != vars[0] || presoldata->setppchashdatas[c].vars[1] != vars[1] )
            {
               removeentry = TRUE;
            }
         }

         if( removeentry )
         {
            /* remove constraint from setppc-hashtable */
            assert(SCIPhashtableExists(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c].cons));
            SCIP_CALL( SCIPhashtableRemove(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c].cons) );

            /* remove hashdata entry from hashtable */
            SCIP_CALL( SCIPhashtableRemove(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[c]) );

            /* release old constraints */
            SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->setppchashdatas[c].cons)) );

            /* release variables */
            SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c].vars[0])) );
            SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c].vars[1])) );

            /* free memory for variables */
            SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[c].vars), 2);

            if( c < presoldata->nsetppchashdatas - 1 )
            {
               /* remove old hashdata entry from hashtable */
               SCIP_CALL( SCIPhashtableRemove(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1]) );
            }

            /* move last content to free position */
            presoldata->setppchashdatas[c].cons = presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1].cons;
            presoldata->setppchashdatas[c].vars = presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1].vars;
            presoldata->setppchashdatas[c].nvars = presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1].nvars;

            if( c < presoldata->nsetppchashdatas - 1 )
            {
               /* add new hashdata entry from hashtable */
               SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[c]) );
            }
            --(presoldata->nsetppchashdatas);
         }
      }

#ifndef NDEBUG
      for( c = presoldata->nsetppchashdatas - 1; c >= 0; --c )
      {
         assert(presoldata->setppchashdatas[c].nvars == 2);
         assert(presoldata->setppchashdatas[c].vars != NULL);
         assert(presoldata->setppchashdatas[c].vars[0] != NULL);
         assert(presoldata->setppchashdatas[c].vars[1] != NULL);
         assert(presoldata->setppchashdatas[c].cons != NULL);
         assert(SCIPconsIsActive(presoldata->setppchashdatas[c].cons));
         assert(SCIPhashtableExists(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[c]));
         assert(SCIPhashtableExists(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c].cons));
      }
#endif
   }

   return SCIP_OKAY;
}

/** refresh useful set-packing information, delete redundant constraints and add new constraints */
static
SCIP_RETCODE correctPresoldata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< data object of presolver */
   SCIP_CONS**           setppcs,            /**< active setppc constraints */
   int                   nsetppcs,           /**< number of active setppc constraints */
   SCIP_CONS**           logicors,           /**< active setppc constraints */
   int                   nlogicors           /**< number of active setppc constraints */
   )
{
   int oldnsetppchashdatas;
   int c;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(setppcs != NULL);
   assert(nsetppcs > 0);
   assert(logicors != NULL);
   assert(nlogicors > 0);
   assert(presoldata->initialized);
   assert(presoldata->setppchashtable != NULL);
   assert(presoldata->logicorhashtable != NULL);

   /* check if there already exist some set-packing and some logicor constraints with the right amount of variables */
   if( !presoldata->usefulsetppcexist || !presoldata->usefullogicorexist )
   {
      SCIP_Bool usefullogicorexisted = presoldata->usefullogicorexist;

      SCIP_CALL( createPresoldata(scip, presoldata, setppcs, nsetppcs, logicors, nlogicors) );

      /* if we already had useful logicor constraints but did not find any useful setppc constraint, the maximal number
       * of variables appearing in a logicor constraint was not updated, so we do it here
       */
      if( usefullogicorexisted && !presoldata->usefulsetppcexist )
      {
         /* correct maximal number of varables in logicor constraints */
         for( c = nlogicors - 1; c >= 0; --c )
         {
            assert(SCIPconsIsActive(logicors[c]));

            /* update maximal entries in a logicor constraint */
            if( presoldata->maxnvarslogicor < SCIPgetNVarsLogicor(scip, logicors[c]) )
               presoldata->maxnvarslogicor = SCIPgetNVarsLogicor(scip, logicors[c]);
         }
      }

      /* no correct logicor or set-packing constraints available, so abort */
      if( !presoldata->usefulsetppcexist || !presoldata->usefullogicorexist )
         return SCIP_OKAY;
   }

   /* correct old data */
   SCIP_CALL( cleanupHashDatas(scip, presoldata) );

   oldnsetppchashdatas = presoldata->nsetppchashdatas;

   /* first update setppc part */
   /* add new setppc constraints */
   for( c = nsetppcs - 1; c >= 0; --c )
   {
      assert(SCIPconsIsActive(setppcs[c]));

      if( SCIPgetTypeSetppc(scip, setppcs[c]) == SCIP_SETPPCTYPE_PACKING && SCIPgetNVarsSetppc(scip, setppcs[c]) == 2 && !SCIPconsIsModifiable(setppcs[c]) )
      {
         /* check if constraint is new, and correct array size if necessary */
         if( !SCIPhashtableExists(presoldata->setppchashtable, (void*) setppcs[c]) )
         {
            SCIP_VAR** setppcvars;
            SCIP_Bool negated[2];

            /* resize array if necessary */
            if( presoldata->nsetppchashdatas == presoldata->ssetppchashdatas )
            {
               int newsize;
               int d;

               newsize = SCIPcalcMemGrowSize(scip, presoldata->nsetppchashdatas + 1);

               /* array already at maximal size */
               if( newsize <= presoldata->ssetppchashdatas )
                  return SCIP_NOMEMORY;

               /* correct hashtable, remove old elements */
               SCIPhashtableRemoveAll(presoldata->hashdatatable);

               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(presoldata->setppchashdatas), presoldata->ssetppchashdatas, newsize) );
               presoldata->ssetppchashdatas = newsize;

               /* add all elements to the hashtable again */
               for( d = presoldata->nsetppchashdatas - 1; d >= 0; --d )
               {
                  SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[d]) );
               }
            }

            /* insert new element in hashtable */
            SCIP_CALL( SCIPhashtableInsert(presoldata->setppchashtable, (void*) setppcs[c]) );

            assert(SCIPgetNVarsSetppc(scip, setppcs[c]) == 2);
            setppcvars = SCIPgetVarsSetppc(scip, setppcs[c]);

            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars), setppcvars, 2) );
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0], &(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0]), &(negated[0])) );
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1], &(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1]), &(negated[1])) );

            if( SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0]) == SCIP_VARSTATUS_MULTAGGR
                  || SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1]) == SCIP_VARSTATUS_MULTAGGR )
            {
               SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars), 2);
               continue;
            }

            presoldata->setppchashdatas[presoldata->nsetppchashdatas].nvars = 2;

            /* capture variables */
            SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0]) );
            SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1]) );

            /* order the variables after their index */
            if( SCIPvarGetIndex(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0]) > SCIPvarGetIndex(presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1]) )
            {
               SCIP_VAR* tmp = presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0];
               presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[0] = presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1];
               presoldata->setppchashdatas[presoldata->nsetppchashdatas].vars[1] = tmp;
            }

            presoldata->setppchashdatas[presoldata->nsetppchashdatas].cons = setppcs[c];

            SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[presoldata->nsetppchashdatas]) );
            SCIP_CALL( SCIPcaptureCons(scip, setppcs[c]) );

            ++(presoldata->nsetppchashdatas);
         }
      }
   }

   /* if we found new set-packing constraints, we want to check against all logicors */
   if( oldnsetppchashdatas < presoldata->nsetppchashdatas )
      presoldata->newsetppchashdatas = TRUE;

   /* now logicor part */
   /* removed last deleted logicor constraints from local presolver data */
   while( presoldata->nusefullogicor > 0 && !SCIPconsIsActive(presoldata->usefullogicor[presoldata->nusefullogicor - 1]) )
   {
      SCIP_CALL( SCIPhashtableRemove(presoldata->logicorhashtable, (void*) presoldata->usefullogicor[presoldata->nusefullogicor - 1]) );
      SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->usefullogicor[presoldata->nusefullogicor - 1])) );

      --(presoldata->nusefullogicor);
   }

   /* remove old inactive logicor constraints */
   for( c = presoldata->nusefullogicor - 1; c >= 0; --c )
   {
      /* update maximal entries in a logicor constraint */
      if( presoldata->maxnvarslogicor < SCIPgetNVarsLogicor(scip, presoldata->usefullogicor[c]) )
         presoldata->maxnvarslogicor = SCIPgetNVarsLogicor(scip, presoldata->usefullogicor[c]);

      if( !SCIPconsIsActive(presoldata->usefullogicor[c]) || SCIPconsIsModifiable(presoldata->usefullogicor[c]) || SCIPgetNVarsLogicor(scip, presoldata->usefullogicor[c]) < 3 )
      {
         SCIP_CALL( SCIPhashtableRemove(presoldata->logicorhashtable, (void*) presoldata->usefullogicor[c]) );
         SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->usefullogicor[c])) );

         presoldata->usefullogicor[c] = presoldata->usefullogicor[presoldata->nusefullogicor - 1];
         --(presoldata->nusefullogicor);
      }
   }

   presoldata->firstchangedlogicor = presoldata->nusefullogicor;
   assert(presoldata->firstchangedlogicor >= 0);

   /* add new logicor constraints */
   for( c = nlogicors - 1; c >= 0; --c )
   {
      assert(SCIPconsIsActive(logicors[c]));

      if( !SCIPconsIsModifiable(logicors[c]) && SCIPgetNVarsLogicor(scip, logicors[c]) >= 3 )
      {
         /* check if constraint is new, and correct array size if necessary */
         if( !SCIPhashtableExists(presoldata->logicorhashtable, (void*) logicors[c]) )
         {
            /* resize array if necessary */
            if( presoldata->nusefullogicor == presoldata->susefullogicor )
            {
               int newsize;

               newsize = SCIPcalcMemGrowSize(scip, presoldata->nusefullogicor + 1);

               /* array already at maximal size */
               if( newsize <= presoldata->susefullogicor )
                  return SCIP_NOMEMORY;

               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(presoldata->usefullogicor), presoldata->susefullogicor, newsize) );
               presoldata->susefullogicor = newsize;
            }

            /* insert new element in hashtable */
            SCIP_CALL( SCIPhashtableInsert(presoldata->logicorhashtable, (void*) logicors[c]) );
            SCIP_CALL( SCIPcaptureCons(scip, logicors[c]) );

            presoldata->usefullogicor[presoldata->nusefullogicor] = logicors[c];
            ++(presoldata->nusefullogicor);

            /* update maximal entries in a logicor constraint */
            if( presoldata->maxnvarslogicor < SCIPgetNVarsLogicor(scip, logicors[c]) )
               presoldata->maxnvarslogicor = SCIPgetNVarsLogicor(scip, logicors[c]);
         }
      }
   }

   return SCIP_OKAY;
}


/** extract and-constraints and set-partitioning constraints */
static
SCIP_RETCODE extractGates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< data object of presolver */
   int                   pos,                /**< position of logicor in usefullogicor array to presolve */
   SCIP_HASHMAP*         varmap,             /**< variable map mapping inactive variables to their active representation */
   SCIP_CONS**           gateconss,          /**< allocated memory for all gate-constraints */
   SCIP_VAR**            activevars,         /**< allocated memory for active variables */
   SCIP_VAR**            posresultants,      /**< allocated memory for all possible resultant variables */
   HASHDATA*             hashdata,           /**< allocated memory for a hashdata object */
   int*                  ndelconss,          /**< pointer to store number of deleted constraints */
   int*                  naddconss           /**< pointer to store number of added constraints */
   )
{
   SCIP_VAR** logicorvars;
   HASHDATA* hashmaphashdata;
   SCIP_CONS* logicor;
   SCIP_Bool negated;
   int ngateconss;
   int nlogicorvars;
   int nposresultants;
   int d;
   int v;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(0 <= pos && pos < presoldata->nusefullogicor);
   assert(gateconss != NULL);
   assert(activevars != NULL);
   assert(posresultants != NULL);
   assert(hashdata != NULL);
   assert(hashdata->vars != NULL);
   assert(hashdata->nvars == 2);
   assert(hashdata->cons == NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);

   assert(presoldata->usefullogicor != NULL);
   logicor = presoldata->usefullogicor[pos];
   assert(logicor != NULL);

   if( !SCIPconsIsActive(logicor) )
      return SCIP_OKAY;

   assert(!SCIPconsIsModifiable(logicor));

   nlogicorvars = SCIPgetNVarsLogicor(scip, logicor);
   assert(nlogicorvars >= 3 && nlogicorvars <= presoldata->maxnvarslogicor);

   logicorvars = SCIPgetVarsLogicor(scip, logicor);
   assert(logicorvars != NULL);

   nposresultants = 0;

   /* get active logicor variables and determine all possible resultants */
   for( d = nlogicorvars - 1; d >= 0; --d )
   {
      /* do not work with fixed variables */
      if( SCIPvarGetLbLocal(logicorvars[d]) > 0.5 || SCIPvarGetUbLocal(logicorvars[d]) < 0.5 )
         return SCIP_OKAY;

      activevars[d] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, logicorvars[d]);

      if( activevars[d] == NULL )
      {
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, logicorvars[d], &(activevars[d]), &negated) );
         SCIP_CALL( SCIPhashmapInsert(varmap, logicorvars[d], activevars[d]) );
      }

      /* determine possible resultants a check if the other variables can appear in a set-packing constraint */
      if( SCIPvarIsNegated(activevars[d]) )
      {
         assert(SCIPvarIsActive(SCIPvarGetNegatedVar(activevars[d])));

         if( SCIPvarGetNLocksDown(SCIPvarGetNegatedVar(activevars[d])) >= nlogicorvars - 1 )
         {
            posresultants[nposresultants] = activevars[d];
            ++nposresultants;
         }
         else if( SCIPvarGetNLocksDown(SCIPvarGetNegatedVar(activevars[d])) == 0 )
            return SCIP_OKAY;
      }
      else
      {
         assert(SCIPvarIsActive(activevars[d]));

         if( SCIPvarGetNLocksUp(activevars[d]) >= nlogicorvars - 1 )
         {
            posresultants[nposresultants] = activevars[d];
            ++nposresultants;
         }
         else if( SCIPvarGetNLocksUp(activevars[d]) == 0 )
            return SCIP_OKAY;
      }
   }

   if( nposresultants == 0 )
      return SCIP_OKAY;

   /* sort variables after indices */
   SCIPsortPtr((void**)activevars, SCIPvarComp, nlogicorvars);

   /* check that we have really different variables, if not remove the constraint from the hashmap and the data
    * storage
    */
   for( d = nlogicorvars - 1; d > 0; --d )
   {
      if( SCIPvarGetIndex(activevars[d]) == SCIPvarGetIndex(activevars[d - 1]) )
      {
         assert(presoldata->usefullogicor[pos] == logicor);

         SCIP_CALL( SCIPhashtableRemove(presoldata->logicorhashtable, (void*) logicor) );
         SCIP_CALL( SCIPreleaseCons(scip, &logicor) );

         presoldata->usefullogicor[pos] = presoldata->usefullogicor[presoldata->nusefullogicor - 1];
         --(presoldata->nusefullogicor);

         return SCIP_OKAY;
      }
   }

   ngateconss = 0;

   for( d = nposresultants - 1; d >= 0; --d )
   {
      ngateconss = 0;

      for( v = nlogicorvars - 1; v >= 0; --v )
      {
         if( activevars[v] == posresultants[d] )
            continue;

         /* variables need to be sorted */
         if( SCIPvarCompare(posresultants[d], activevars[v]) > 0 )
         {
            hashdata->vars[0] = activevars[v];
            hashdata->vars[1] = posresultants[d];
         }
         else
         {
            hashdata->vars[0] = posresultants[d];
            hashdata->vars[1] = activevars[v];
         }

         hashmaphashdata = (HASHDATA*) SCIPhashtableRetrieve(presoldata->hashdatatable, (void*) hashdata);

         if( hashmaphashdata != NULL && SCIPconsIsActive(hashmaphashdata->cons) )
         {
            gateconss[ngateconss] = hashmaphashdata->cons;
            ++ngateconss;
         }
         else
            break;
      }
      if( ngateconss == nlogicorvars - 1 )
         break;
   }

   /* @todo, check for clique of all variables except the resultant */
   /* check if we have a set-partitioning 'gate' */
   if( ngateconss == nlogicorvars - 1 && nlogicorvars == 3 )
   {
      assert(d >= 0 && d < nposresultants);
      assert(ngateconss >= 2);

      if( activevars[0] == posresultants[d] )
      {
         hashdata->vars[0] = activevars[1];
         hashdata->vars[1] = activevars[2];
      }
      else if( activevars[1] == posresultants[d] )
      {
         hashdata->vars[0] = activevars[0];
         hashdata->vars[1] = activevars[2];
      }
      else
      {
         assert(activevars[2] == posresultants[d]);
         hashdata->vars[0] = activevars[0];
         hashdata->vars[1] = activevars[1];
      }

      hashmaphashdata = (HASHDATA*) SCIPhashtableRetrieve(presoldata->hashdatatable, (void*) hashdata);
      assert(hashmaphashdata == NULL || hashmaphashdata->cons != NULL);

      if( hashmaphashdata != NULL && SCIPconsIsActive(hashmaphashdata->cons) )
      {
         gateconss[ngateconss] = hashmaphashdata->cons;
         ++ngateconss;
      }
   }

   /* did we find enough (>= number of variables in logicor - 1) set-packing constraints for an upgrade to either
    * an and-constraint or even a set-partitioning constraint
    */
   if( ngateconss == nlogicorvars || (ngateconss >= nlogicorvars - 1 && !presoldata->onlysetpart))
   {
      SCIP_CONS* newcons;
      char name[SCIP_MAXSTRLEN];
      SCIP_Bool initial;
      SCIP_Bool separate;
      SCIP_Bool enforce;
      SCIP_Bool check;
      SCIP_Bool propagate;
      SCIP_Bool local;
      SCIP_Bool modifiable;
      SCIP_Bool dynamic;
      SCIP_Bool removable;
      SCIP_Bool stickingatnode;
      int i;

      assert(ngateconss <= nlogicorvars);
      assert(d >= 0 && d < nposresultants);

      initial = SCIPconsIsInitial(logicor);
      separate = SCIPconsIsSeparated(logicor);
      enforce = SCIPconsIsEnforced(logicor);
      check = SCIPconsIsChecked(logicor);
      propagate = SCIPconsIsPropagated(logicor);
      local = SCIPconsIsLocal(logicor);
      modifiable = SCIPconsIsModifiable(logicor);
      dynamic = SCIPconsIsDynamic(logicor);
      removable = SCIPconsIsRemovable(logicor);
      stickingatnode = SCIPconsIsStickingAtNode(logicor);

#ifdef SCIP_DEBUG
      if( ngateconss == nlogicorvars )
      {
         SCIPdebugMsg(scip, "Following constraints form a set-partitioning constraint.\n");
      }
      else
      {
         SCIPdebugMsg(scip, "Following constraints form an and-constraint.\n");
      }
#endif

      for( v = ngateconss - 1; v >= 0; --v )
      {
         assert(gateconss[v] != NULL);

         initial |= SCIPconsIsInitial(gateconss[v]);
         separate |= SCIPconsIsSeparated(gateconss[v]);
         enforce |= SCIPconsIsEnforced(gateconss[v]);
         check |= SCIPconsIsChecked(gateconss[v]);
         propagate |= SCIPconsIsPropagated(gateconss[v]);
         local &= SCIPconsIsLocal(gateconss[v]);
         modifiable &= SCIPconsIsModifiable(gateconss[v]);
         dynamic &= SCIPconsIsDynamic(gateconss[v]);
         removable &= SCIPconsIsRemovable(gateconss[v]);
         stickingatnode &= SCIPconsIsStickingAtNode(gateconss[v]);

         SCIPdebugPrintCons(scip, gateconss[v], NULL);

         SCIP_CALL( SCIPdelCons(scip, gateconss[v]) );
         ++(*ndelconss);
      }

      SCIPdebugPrintCons(scip, logicor, NULL);

      if( ngateconss == nlogicorvars - 1 )
      {
         SCIP_VAR** consvars;

         assert(!presoldata->onlysetpart);

         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, ngateconss) );
         i = 0;

         /* determine and operands */
         for( v = nlogicorvars - 1; v >= 0; --v )
         {
            if( activevars[v] == posresultants[d] )
               continue;

            SCIP_CALL( SCIPgetNegatedVar(scip, activevars[v], &consvars[i]) );
            ++i;
         }
         assert(i == ngateconss);

         /* create and add "and" constraint for the extracted gate */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "andgate_%d", presoldata->ngates);
         SCIP_CALL( SCIPcreateConsAnd(scip, &newcons, name, posresultants[d], ngateconss, consvars,
               initial, separate, enforce, check, propagate,
               local, modifiable, dynamic, removable, stickingatnode) );

         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugMsg(scip, "-------------->\n");
         SCIPdebugPrintCons(scip, newcons, NULL);

         ++(*naddconss);
         ++(presoldata->ngates);

         SCIP_CALL( SCIPdelCons(scip, logicor) );
         ++(*ndelconss);

         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

         SCIPfreeBufferArray(scip, &consvars);
      }
      else
      {
         assert(ngateconss == nlogicorvars);

         /* create and add set-partitioning constraint */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "setpart_%d", presoldata->ngates);
         SCIP_CALL( SCIPcreateConsSetpart(scip, &newcons, name, nlogicorvars, activevars,
               initial, separate, enforce, check, propagate,
               local, modifiable, dynamic, removable, stickingatnode) );

         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIPdebugMsg(scip, "-------------->\n");
         SCIPdebugPrintCons(scip, newcons, NULL);

         ++(*naddconss);
         ++(presoldata->ngates);

         SCIP_CALL( SCIPdelCons(scip, logicor) );
         ++(*ndelconss);

         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyGateextraction)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolGateextraction(scip) );

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeGateextraction)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   if( presoldata->hashdatatable != NULL )
   {
      assert(presoldata->setppchashtable != NULL);
      assert(presoldata->logicorhashtable != NULL);

      SCIPhashtableFree(&(presoldata->logicorhashtable));
      SCIPhashtableFree(&(presoldata->setppchashtable));
      SCIPhashtableFree(&(presoldata->hashdatatable));
   }

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitGateextraction)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   int c;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* release old constraints */
   for( c = presoldata->nusefullogicor - 1; c >= 0; --c )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->usefullogicor[c])) );
   }

   if( presoldata->usefullogicorexist )
   {
      SCIPfreeBlockMemoryArray(scip, &presoldata->usefullogicor, presoldata->susefullogicor);
   }

   if( presoldata->usefulsetppcexist )
   {
      assert(presoldata->setppchashdatas != NULL || presoldata->nsetppchashdatas == 0);
      for( c = presoldata->nsetppchashdatas - 1; c >= 0; --c )
      {
         assert(presoldata->setppchashdatas[c].cons != NULL);
         assert(presoldata->setppchashdatas[c].vars != NULL);

         /* remove constraint from setppc-hashtable */
         assert(SCIPhashtableExists(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c].cons));
         SCIP_CALL( SCIPhashtableRemove(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c].cons) );


         /* remove hashdata entry from hashtable */
         SCIP_CALL( SCIPhashtableRemove(presoldata->hashdatatable, (void*) &presoldata->setppchashdatas[c]) );

         /* release old constraints */
         SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->setppchashdatas[c].cons)) );

         /* release variables */
         SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c].vars[0])) );
         SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c].vars[1])) );

         /* free memory for variables */
         SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[c].vars), 2);
      }

      SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas), presoldata->ssetppchashdatas);
   }

   if( presoldata->hashdatatable != NULL )
   {
      assert(presoldata->setppchashtable != NULL);
      assert(presoldata->logicorhashtable != NULL);

      /* clear old hashtable entries */
      SCIPhashtableRemoveAll(presoldata->hashdatatable);
      SCIPhashtableRemoveAll(presoldata->setppchashtable);
      SCIPhashtableRemoveAll(presoldata->logicorhashtable);
   }

   presoldata->nusefullogicor = 0;
   presoldata->susefullogicor = 0;
   presoldata->nsetppchashdatas = 0;
   presoldata->ssetppchashdatas = 0;
   presoldata->firstchangedlogicor = -1;
   presoldata->ngates = 0;
   presoldata->usefullogicorexist = FALSE;
   presoldata->usefulsetppcexist = FALSE;
   presoldata->newsetppchashdatas = FALSE;
   presoldata->initialized = FALSE;

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreGateextraction)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreGateextraction)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecGateextraction)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_HASHMAP* varmap;
   HASHDATA  hashdata;
   SCIP_VAR* tmpvars[2];
   SCIP_CONSHDLR* conshdlrsetppc;
   SCIP_CONSHDLR* conshdlrlogicor;
   SCIP_CONSHDLR* conshdlrand;
   SCIP_CONS** setppcconss;
   SCIP_CONS** logicorconss;
   int nsetppcconss;
   int nlogicorconss;
   int size;
   int c;
   SCIP_Bool paramvalue;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

#if 0 /* need to include cons_knapsack on top of this file */
   /* check for possible knapsacks that form with a logicor a weak relaxation of an and-constraint
    *
    * the weak relaxation of an and-constraint looks like:
    *   - row1:             resvar - v1 - ... - vn >= 1-n
    *   - row2:           n*resvar - v1 - ... - vn <= 0.0
    *
    * which look like the following contraints
    *   - logicor:          resvar + ~v1 + ... + ~vn >= 1
    *   - knapsack:       n*resvar + ~v1 + ... + ~vn <= n
    */
   {
      SCIP_CONSHDLR* conshdlrknapsack;
      SCIP_CONS** knapsackconss;
      int nknapsackconss;
      SCIP_VAR** vars;
      SCIP_Longint* vals;
      SCIP_Longint capacity;
      int nvars;

      conshdlrknapsack = SCIPfindConshdlr(scip, "knapsack");

      /* get number of active constraints */
      knapsackconss = SCIPconshdlrGetConss(conshdlrknapsack);
      nknapsackconss = SCIPconshdlrGetNActiveConss(conshdlrknapsack);
      assert(nknapsackconss >= 0);
      assert(knapsackconss != NULL || nknapsackconss == 0);

      for( c = nknapsackconss - 1; c >= 0; --c )
      {
         /* not implemented in master branch, but the constraint may be already sorted */
         /*SCIPsortKnapsack(scip, knapsackconss[c]);*/

         nvars = SCIPgetNVarsKnapsack(scip, knapsackconss[c]);
         vals = SCIPgetWeightsKnapsack(scip, knapsackconss[c]);
         vars = SCIPgetVarsKnapsack(scip, knapsackconss[c]);
         capacity = SCIPgetCapacityKnapsack(scip, knapsackconss[c]);

         if( nvars > 1 && capacity == nvars - 1 && vals[0] == capacity && vals[1] == 1 )
         {
            printf("possible knapsack for gate extraction\n");
         }
      }
   }
#endif

   /* get necessary constraint handlers */
   conshdlrsetppc = SCIPfindConshdlr(scip, "setppc");
   conshdlrlogicor = SCIPfindConshdlr(scip, "logicor");

   if( conshdlrsetppc == NULL || conshdlrlogicor == NULL )
      return SCIP_OKAY;

   /* get number of active constraints */
   nsetppcconss = SCIPconshdlrGetNActiveConss(conshdlrsetppc);
   assert(nsetppcconss >= 0);
   nlogicorconss = SCIPconshdlrGetNActiveConss(conshdlrlogicor);
   assert(nlogicorconss >= 0);

   if( nsetppcconss == 0 || nlogicorconss == 0 )
      return SCIP_OKAY;

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   conshdlrand = SCIPfindConshdlr(scip, "and");

   /* need and-constraint handler to extract and-gates */
   if( conshdlrand == NULL )
   {
      /* nothing to do when we cannot extract anything */
      if( !presoldata->searchequations )
         return SCIP_OKAY;
      else
      {
         /* make sure that we correct the parameter for only extrating set-partitioning constraints */
         if( SCIPisParamFixed(scip, "presolving/" PRESOL_NAME "/onlysetpart") )
         {
            SCIPwarningMessage(scip, "unfixing parameter <presolving/" PRESOL_NAME "/onlysetpart> in gate extration presolver\n");
            SCIP_CALL( SCIPunfixParam(scip, "presolving/" PRESOL_NAME "/onlysetpart") );
         }
         SCIP_CALL( SCIPsetBoolParam(scip, "presolving/" PRESOL_NAME "/onlysetpart", TRUE) );
         assert(presoldata->onlysetpart);
      }
   }

   paramvalue = FALSE;
   if( conshdlrand != NULL && SCIPgetBoolParam(scip, "constraints/and/linearize", &paramvalue) == SCIP_OKAY )
   {
      if( paramvalue )
      {
	 SCIPwarningMessage(scip, "Gate-presolving is the 'counterpart' of linearizing all and-constraints, so enabling both presolving steps at ones does not make sense.\n");
      }
   }
   *result = SCIP_DIDNOTFIND;

   /* get active constraints */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &setppcconss, SCIPconshdlrGetConss(conshdlrsetppc), nsetppcconss) ); /*lint !e666*/

   assert(setppcconss != NULL);
   logicorconss = SCIPconshdlrGetConss(conshdlrlogicor);
   assert(logicorconss != NULL);

   /* first we need to initialized the hashtables if not yet done */
   if( presoldata->hashdatatable == NULL )
   {
      SCIP_CALL( presoldataInitHashtables(scip, presoldata) );
   }
   assert(presoldata->hashdatatable != NULL);
   assert(presoldata->setppchashtable != NULL);
   assert(presoldata->logicorhashtable != NULL);

   presoldata->newsetppchashdatas = FALSE;

   if( !presoldata->initialized )
   {
      assert(presoldata->usefullogicor == NULL);

      /* create useful set-packing information by adding new set-packing constraints with two variables */
      SCIP_CALL( createPresoldata(scip, presoldata, setppcconss, nsetppcconss, logicorconss, nlogicorconss) );
   }
   else
   {
      /* refresh useful set-packing information, delete redundant constraints and add new constraints */
      SCIP_CALL( correctPresoldata(scip, presoldata, setppcconss, nsetppcconss, logicorconss, nlogicorconss) );
   }
   assert(presoldata->initialized);

   if( presoldata->nusefullogicor == 0 )
      goto TERMINATE;

   /* move the biggate extraction to front or back by sort the logicors after number of variables */

   if( presoldata->sorting != 0 )
   {
      int* lengths;

      SCIP_CALL( SCIPallocBufferArray(scip, &lengths, presoldata->nusefullogicor) );

      for( c = presoldata->nusefullogicor - 1; c >= 0; --c )
      {
         lengths[c] = SCIPgetNVarsLogicor(scip, presoldata->usefullogicor[c]);
      }

      if( presoldata->sorting == -1 )
         SCIPsortDownIntPtr(lengths, (void**)presoldata->usefullogicor, presoldata->nusefullogicor);
      else
         SCIPsortIntPtr(lengths, (void**)presoldata->usefullogicor, presoldata->nusefullogicor);

      SCIPfreeBufferArray(scip, &lengths);
   }

   /* maximal number of binary variables */
   size = SCIPgetNBinVars(scip) + SCIPgetNImplVars(scip);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), size) );

   /* search for set-partitioning constraints arising from a logicor and a set-packing constraints with equal variables */
   if( presoldata->searchequations && !SCIPisStopped(scip) )
   {
      SCIP_HASHTABLE* setppchashdatatable;
      HASHDATA** setppchashdatas;
      HASHDATA* setppchashdatastore;
      HASHDATA* hashmaphashdata;
      SCIP_CONS* logicor;
      SCIP_CONS* setppc;
      SCIP_VAR** logicorvars;
      SCIP_VAR** setppcvars;
      SCIP_VAR** activevarslogicor;
      SCIP_VAR** activevarssetppc;
      SCIP_Bool negated;
      int nsetppchashdatas;
      int nlogicorvars;
      int nsetppcvars;
      int d;
      int v;

      assert(nsetppcconss > 0);

      /* create local hashtable */
      SCIP_CALL( SCIPhashtableCreate(&setppchashdatatable, SCIPblkmem(scip), nsetppcconss, SCIPhashGetKeyStandard, setppcHashdataKeyEqCons, setppcHashdataKeyValCons, (void*) scip) );

      /* maximal number of binary variables */
      size = presoldata->maxnvarslogicor;
      assert(size >= 3);

      /* get temporary memory */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &setppchashdatastore, nsetppcconss) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &setppchashdatas, nsetppcconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &activevarssetppc, size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &activevarslogicor, size) );

      hashdata.cons = NULL;

      nsetppchashdatas = 0;

      /* collect all set-packing/-partitioning constraints and corresponding data to be able to search faster */
      for( d = nsetppcconss - 1; d >= 0; --d )
      {
         setppc = setppcconss[d];
         assert(setppc != NULL);

         if( SCIPconsIsDeleted(setppc) )
            continue;

         /* @todo if of interest could also be implemented for set-covering constraints */
#if 1
         if( SCIPgetTypeSetppc(scip, setppc) == SCIP_SETPPCTYPE_COVERING )
            continue;
#endif

         nsetppcvars = SCIPgetNVarsSetppc(scip, setppc);

         if( nsetppcvars < 2 )
            continue;

         if( SCIPconsIsModifiable(setppc) )
            continue;

         /* to big setppc constraints are picked out */
         if( nsetppcvars > size )
            continue;

         setppcvars = SCIPgetVarsSetppc(scip, setppc);
         assert(setppcvars != NULL);

         /* get active setppc variables */
         for( v = nsetppcvars - 1; v >= 0; --v )
         {
            /* do not work with fixed variables */
            if( SCIPvarGetLbLocal(setppcvars[v]) > 0.5 || SCIPvarGetUbLocal(setppcvars[v]) < 0.5 )
               break;

            activevarssetppc[v] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, setppcvars[v]);

            if( activevarssetppc[v] == NULL )
            {
               SCIP_CALL( SCIPgetBinvarRepresentative(scip, setppcvars[v], &(activevarssetppc[v]), &negated) );
               SCIP_CALL( SCIPhashmapInsert(varmap, setppcvars[v], activevarssetppc[v]) );
            }
         }

         /* if we found a fixed variable we want disregard this constraint */
         if( v >= 0 )
            continue;

         /* variables need to be sorted after indices to be able to do a fast comparison */
         SCIPsortPtr((void**)activevarssetppc, SCIPvarComp, nsetppcvars);

         setppchashdatas[nsetppchashdatas] = &(setppchashdatastore[nsetppchashdatas]);

         /* memorize set-packing data */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(setppchashdatas[nsetppchashdatas]->vars), activevarssetppc, nsetppcvars) );

         setppchashdatas[nsetppchashdatas]->nvars = nsetppcvars;
         setppchashdatas[nsetppchashdatas]->cons = setppc;
         /* need to capture this constraint, because it might get deleted during the process */
         SCIP_CALL( SCIPcaptureCons(scip, setppc) );

         /* add entry to local hashtable */
         SCIP_CALL( SCIPhashtableInsert(setppchashdatatable, (void*) setppchashdatas[nsetppchashdatas]) );
         ++nsetppchashdatas;
      }

      /* check all (new) logicors against all collected set-packing/-partitioning constraints */
      for( c = nlogicorconss - 1; c >= 0 && !SCIPisStopped(scip); --c )
      {
         logicor = logicorconss[c];
         assert(logicor != NULL);

         if( SCIPconsIsDeleted(logicor) )
            continue;

         nlogicorvars = SCIPgetNVarsLogicor(scip, logicor);

         if( nlogicorvars < 2 )
            continue;

         if( SCIPconsIsModifiable(logicor) )
            continue;

         assert(nlogicorvars <= size);

         logicorvars = SCIPgetVarsLogicor(scip, logicor);
         assert(logicorvars != NULL);

         /* get active logicor variables */
         for( v = nlogicorvars - 1; v >= 0; --v )
         {
            /* do not work with fixed variables */
            if( SCIPvarGetLbLocal(logicorvars[v]) > 0.5 || SCIPvarGetUbLocal(logicorvars[v]) < 0.5 )
               break;

            activevarslogicor[v] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, logicorvars[v]);

            /* if image does not exist, then there is no corresponding set-packing constraint */
            if( activevarslogicor[v] == NULL )
               break;
         }

         if( v == -1 )
         {
            /* need sorting to be able to find the correct hashdata element */
            SCIPsortPtr((void**)activevarslogicor, SCIPvarComp, nlogicorvars);

            hashdata.nvars = nlogicorvars;
            hashdata.vars = activevarslogicor;

            hashmaphashdata = (HASHDATA*) SCIPhashtableRetrieve(setppchashdatatable, (void*) &hashdata);
            assert(hashmaphashdata == NULL || hashmaphashdata->cons != NULL);

            if( hashmaphashdata != NULL && !SCIPconsIsDeleted(hashmaphashdata->cons) )
            {
               SCIP_Bool initial;
               SCIP_Bool separate;
               SCIP_Bool enforce;
               SCIP_Bool check;
               SCIP_Bool propagate;
               SCIP_Bool local;
               SCIP_Bool modifiable;
               SCIP_Bool dynamic;
               SCIP_Bool removable;
               SCIP_Bool stickingatnode;

               setppc = hashmaphashdata->cons;
               assert(SCIPconsGetHdlr(setppc) == SCIPfindConshdlr(scip, "setppc"));

               initial = SCIPconsIsInitial(logicor) || SCIPconsIsInitial(setppc);
               separate = SCIPconsIsSeparated(logicor) || SCIPconsIsSeparated(setppc);
               enforce = SCIPconsIsEnforced(logicor) || SCIPconsIsEnforced(setppc);
               check = SCIPconsIsChecked(logicor) || SCIPconsIsChecked(setppc);
               propagate = SCIPconsIsPropagated(logicor) || SCIPconsIsPropagated(setppc);
               local = SCIPconsIsLocal(logicor) && SCIPconsIsLocal(setppc);
               modifiable = SCIPconsIsModifiable(logicor) && SCIPconsIsModifiable(setppc);
               dynamic = SCIPconsIsDynamic(logicor) && SCIPconsIsDynamic(setppc);
               removable = SCIPconsIsRemovable(logicor) && SCIPconsIsRemovable(setppc);
               stickingatnode = SCIPconsIsStickingAtNode(logicor) && SCIPconsIsStickingAtNode(setppc);

               /* check if logicor is redundant against a set-partitioning constraint */
               if( SCIPgetTypeSetppc(scip, setppc) == SCIP_SETPPCTYPE_PARTITIONING )
               {
                  SCIP_CALL( SCIPsetConsInitial(scip, setppc, initial) );
                  SCIP_CALL( SCIPsetConsSeparated(scip, setppc, separate) );
                  SCIP_CALL( SCIPsetConsEnforced(scip, setppc, enforce) );
                  SCIP_CALL( SCIPsetConsChecked(scip, setppc, check) );
                  SCIP_CALL( SCIPsetConsPropagated(scip, setppc, propagate) );
                  SCIP_CALL( SCIPsetConsLocal(scip, setppc, local) );
                  SCIP_CALL( SCIPsetConsModifiable(scip, setppc, modifiable) );
                  SCIP_CALL( SCIPsetConsDynamic(scip, setppc, dynamic) );
                  SCIP_CALL( SCIPsetConsRemovable(scip, setppc, removable) );
                  SCIP_CALL( SCIPsetConsStickingAtNode(scip, setppc, stickingatnode) );

                  SCIPdebugMsg(scip, "Following logicor is redundant to the set-partitioning constraint.\n");
                  SCIPdebugPrintCons(scip, logicor, NULL);
                  SCIPdebugPrintCons(scip, setppc, NULL);
               }
               else
               {
                  SCIP_CONS* newcons;
                  char name[SCIP_MAXSTRLEN];

                  assert(SCIPgetTypeSetppc(scip, setppc) == SCIP_SETPPCTYPE_PACKING);

                  SCIPdebugMsg(scip, "Following logicor and set-packing constraints form a set-partitioning constraint.\n");
                  SCIPdebugPrintCons(scip, logicor, NULL);
                  SCIPdebugPrintCons(scip, setppc, NULL);

                  /* create and add set-partitioning constraint */
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "setpart_%d", presoldata->ngates);
                  SCIP_CALL( SCIPcreateConsSetpart(scip, &newcons, name, nlogicorvars, activevarslogicor,
                        initial, separate, enforce, check, propagate,
                        local, modifiable, dynamic, removable, stickingatnode) );

                  SCIP_CALL( SCIPaddCons(scip, newcons) );
                  SCIPdebugMsg(scip, "-------------->\n");
                  SCIPdebugPrintCons(scip, newcons, NULL);

                  ++(*naddconss);
                  ++(presoldata->ngates);

                  /* delete redundant set-packing constraint */
                  SCIP_CALL( SCIPdelCons(scip, setppc) );
                  ++(*ndelconss);

                  SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
               }

               /* delete redundant logicor constraint */
               SCIP_CALL( SCIPdelCons(scip, logicor) );
               ++(*ndelconss);
            }
         }
      }

      /* need to clear/release parts of hashdata objects */
      for( d = nsetppchashdatas - 1; d >= 0; --d )
      {
         /* need to release captured constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &(setppchashdatas[d]->cons)) );
         /* need to free copied memory */
         SCIPfreeBlockMemoryArray(scip, &(setppchashdatas[d]->vars), setppchashdatas[d]->nvars);
      }

      /* delete local hashtable */
      SCIPhashtableFree(&setppchashdatatable);

      /* free all temporary memory */
      SCIPfreeBufferArray(scip, &activevarslogicor);
      SCIPfreeBufferArray(scip, &activevarssetppc);
      SCIPfreeBlockMemoryArray(scip, &setppchashdatas, nsetppcconss);
      SCIPfreeBlockMemoryArray(scip, &setppchashdatastore, nsetppcconss);
   }

   /* we do not have any useful set-packing or logicor constraint, or since last run did not get any new constraints, so abort */
   if( presoldata->nsetppchashdatas == 0 || (presoldata->firstchangedlogicor == presoldata->nusefullogicor && !presoldata->newsetppchashdatas) )
   {
      SCIPhashmapFree(&varmap);
      goto TERMINATE;
   }

   assert(presoldata->usefullogicor != NULL);
   assert(presoldata->nusefullogicor > 0);
   assert(presoldata->firstchangedlogicor >= 0);
   assert(presoldata->nsetppchashdatas > 0);

   /* search for gates */
   if( presoldata->nsetppchashdatas > 0 && !SCIPisStopped(scip) )
   {
      SCIP_CONS** gateconss;
      SCIP_VAR** activevars;
      SCIP_VAR** posresultants;
      int endloop;

      /* if we found new setppcs we want to check all logicors again */
      if( presoldata->newsetppchashdatas )
	 endloop = 0;
      else
	 endloop = MAX(presoldata->firstchangedlogicor, 0);

      assert(presoldata->maxnvarslogicor >= 3);
      SCIP_CALL( SCIPallocBufferArray(scip, &gateconss, presoldata->maxnvarslogicor) );
      SCIP_CALL( SCIPallocBufferArray(scip, &activevars, presoldata->maxnvarslogicor) );
      SCIP_CALL( SCIPallocBufferArray(scip, &posresultants, presoldata->maxnvarslogicor) );

      hashdata.nvars = 2;
      hashdata.cons = NULL;
      /* assign array of two variables as temporary storage to hashdata */
      hashdata.vars = tmpvars;

      /* check all (new) logicors against all set-packing constraints, to extract and-constraints with two or more
       * operands or set-partitioning constraints three or more variables
       */
      for( c = presoldata->nusefullogicor - 1; c >= endloop && !SCIPisStopped(scip); --c )
      {
         assert(presoldata->usefullogicor[c] != NULL);

         /* logicor constraint has the form: x + y + z >= 1
          *
          * find set-packing constraints:  (~x + ~y >= 1 and ~x + ~z >= 1)  <=>  (x + y <= 1 and x + z <= 1)
          *
          * - these three constraints are aquivalent to: x = ~y * ~z (x = AND(~y,~z))
          *
          * if an additional set-packing constraint exists: y + z <= 1
          *
          * - these four constraints are aquivalent to: x + y + z = 1
          */
         SCIP_CALL( extractGates(scip, presoldata, c, varmap, gateconss, activevars, posresultants, &hashdata, ndelconss, naddconss) );
      }

      SCIPfreeBufferArray(scip, &posresultants);
      SCIPfreeBufferArray(scip, &activevars);
      SCIPfreeBufferArray(scip, &gateconss);
   }

   SCIPhashmapFree(&varmap);

 TERMINATE:
   SCIPfreeBufferArray(scip, &setppcconss);

   /* remove old setppchashdatas objects */
   SCIP_CALL( cleanupHashDatas(scip, presoldata) );

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the gateextraction presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolGateextraction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* alloc presolve data object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* initialize gateextraction presolver data */
   presoldataInit(presoldata);

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecGateextraction, presoldata) );

   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyGateextraction) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeGateextraction) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitGateextraction) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreGateextraction) );
   SCIP_CALL( SCIPsetPresolExitpre(scip, presol, presolExitpreGateextraction) );

   /* add gateextraction presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/onlysetpart",
         "should we only try to extract set-partitioning constraints and no and-constraints",
         &presoldata->onlysetpart, TRUE, DEFAULT_ONLYSETPART, NULL, NULL) );

   /* add gateextraction presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/searchequations",
         "should we try to extract set-partitioning constraint out of one logicor and one corresponding set-packing constraint",
         &presoldata->searchequations, TRUE, DEFAULT_SEARCHEQUATIONS, NULL, NULL) );

   /* add gateextraction presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/sorting",
         "order logicor contraints to extract big-gates before smaller ones (-1), do not order them (0) or order them to extract smaller gates at first (1)",
         &presoldata->sorting, TRUE, DEFAULT_SORTING, -1, 1, NULL, NULL) );

   return SCIP_OKAY;
}
