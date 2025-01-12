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

/**@file   scip_sepa.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for separator plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/scip_sepa.h"
#include "scip/sepa.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/tree.h"

/** creates a separator and includes it in SCIP.
 *
 *  @note method has all separator callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludeSepaBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeSepa(
   SCIP*                 scip,               /**< SCIP data structure */
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
   SCIP_SEPA* sepa;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeSepa", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether separator is already present */
   if( SCIPfindSepa(scip, name) != NULL )
   {
      SCIPerrorMessage("separator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPsepaCreate(&sepa, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, freq, maxbounddist, usessubscip, delay,
         sepacopy, sepafree, sepainit, sepaexit, sepainitsol, sepaexitsol, sepaexeclp, sepaexecsol, sepadata) );
   SCIP_CALL( SCIPsetIncludeSepa(scip->set, sepa) );

   return SCIP_OKAY;
}

/** creates a separator and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetSepaInit(), SCIPsetSepaFree(),
 *  SCIPsetSepaInitsol(), SCIPsetSepaExitsol(), SCIPsetSepaCopy(), SCIPsetExit().
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeSepa() instead
 */
SCIP_RETCODE SCIPincludeSepaBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA**           sepa,               /**< reference to a separator, or NULL */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp)),    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol)),   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   SCIP_SEPA* sepaptr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeSepaBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether separator is already present */
   if( SCIPfindSepa(scip, name) != NULL )
   {
      SCIPerrorMessage("separator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPsepaCreate(&sepaptr, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, freq, maxbounddist, usessubscip, delay,
         NULL, NULL, NULL, NULL, NULL, NULL, sepaexeclp, sepaexecsol, sepadata) );

   assert(sepaptr != NULL);

   SCIP_CALL( SCIPsetIncludeSepa(scip->set, sepaptr) );

   if( sepa != NULL)
      *sepa = sepaptr;

   return SCIP_OKAY;
}

/** sets copy method of separator */
SCIP_RETCODE SCIPsetSepaCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPACOPY    ((*sepacopy))       /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSepaCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(sepa != NULL);

   SCIPsepaSetCopy(sepa, sepacopy);

   return SCIP_OKAY;
}

/** sets destructor method of separator */
SCIP_RETCODE SCIPsetSepaFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAFREE    ((*sepafree))       /**< destructor of separator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSepaFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(sepa != NULL);

   SCIPsepaSetFree(sepa, sepafree);

   return SCIP_OKAY;
}

/** sets initialization method of separator */
SCIP_RETCODE SCIPsetSepaInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINIT    ((*sepainit))       /**< initialize separator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSepaInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(sepa != NULL);

   SCIPsepaSetInit(sepa, sepainit);

   return SCIP_OKAY;
}

/** sets deinitialization method of separator */
SCIP_RETCODE SCIPsetSepaExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit))       /**< deinitialize separator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSepaExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(sepa != NULL);

   SCIPsepaSetExit(sepa, sepaexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of separator */
SCIP_RETCODE SCIPsetSepaInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol))    /**< solving process initialization method of separator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSepaInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(sepa != NULL);

    SCIPsepaSetInitsol(sepa, sepainitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of separator */
SCIP_RETCODE SCIPsetSepaExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol))    /**< solving process deinitialization method of separator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSepaExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(sepa != NULL);

   SCIPsepaSetExitsol(sepa, sepaexitsol);

   return SCIP_OKAY;
}

/** returns the separator of the given name, or NULL if not existing */
SCIP_SEPA* SCIPfindSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of separator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindSepa(scip->set, name);
}

/** returns the array of currently available separators */
SCIP_SEPA** SCIPgetSepas(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortSepas(scip->set);

   return scip->set->sepas;
}

/** returns the number of currently available separators */
int SCIPgetNSepas(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nsepas;
}

/** sets the priority of a separator */
SCIP_RETCODE SCIPsetSepaPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   int                   priority            /**< new priority of the separator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsepaSetPriority(sepa, scip->set, priority);

   return SCIP_OKAY;
}

/** declares separator to be a parent separator
 *
 *  Parent separators generate cuts of several types. To distinguish these cuts, they create child separators, which are
 *  only needed to detect which cuts are applied.
 */
void SCIPsetSepaIsParentsepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(scip != NULL);
   assert(sepa != NULL);

   SCIPsepaSetIsParentsepa(sepa);
}

/** sets the parent separator
 *
 *  Informs SCIP that the separator @p sepa depends on the parent separator @p parentsepa.
 */
void SCIPsetSepaParentsepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPA*            parentsepa          /**< parent separator */
   )
{
   assert(scip != NULL);
   assert(sepa != NULL);

   SCIPsepaSetParentsepa(sepa, parentsepa);
}

#undef SCIPgetSepaMinEfficacy

/** gets value of minimal efficacy for a cut to enter the LP
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return value of "separating/minefficacyroot" if at root node, otherwise value of "separating/minefficacy"
 */
SCIP_Real SCIPgetSepaMinEfficacy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->tree != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSepaMinEfficacy", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeGetCurrentDepth(scip->tree) != 0 )
      return scip->set->sepa_minefficacyroot;
   return scip->set->sepa_minefficacy;
}
