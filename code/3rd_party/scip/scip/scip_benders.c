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

/**@file   scip_benders.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for Benders decomposition
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

#include "scip/benderscut.h"
#include "scip/benders.h"
#include "scip/cons_linear.h"
#include "scip/debug.h"
#include "scip/dcmp.h"
#include "scip/pub_benders.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/scip.h"
#include "scip/scip_message.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates a Benders' decomposition and includes it in SCIP
 *
 *  To use the Benders' decomposition for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note method has all Benders' decomposition callbacks as arguments and is thus changed every time a new callback is
 *        added in future releases; consider using SCIPincludeBendersBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve)),/**< the execution method of the Benders' decomposition algorithm */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve)),/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub)),/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   SCIP_BENDERS* benders;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBenders", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether pricer is already present */
   if( SCIPfindBenders(scip, name) != NULL )
   {
      SCIPerrorMessage("benders <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benderssolvesubconvex == NULL && benderssolvesub == NULL && bendersfreesub != NULL)
      || ((benderssolvesubconvex != NULL || benderssolvesub != NULL) && bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that if bendersFreesub%s is "
         "implemented at least one of bendersSolvesubconvex%s or bendersSolvesub%s are implemented, "
         "or if bendersFreesub%s is not implemented, then none are implented.\n", name, name, name, name, name);
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPbendersCreate(&benders, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         cutlp, cutpseudo, cutrelax, shareauxvars, benderscopy, bendersfree, bendersinit, bendersexit, bendersinitpre,
         bendersexitpre, bendersinitsol, bendersexitsol, bendersgetvar, benderscreatesub, benderspresubsolve,
         benderssolvesubconvex, benderssolvesub, benderspostsolve, bendersfreesub, bendersdata) );
   SCIP_CALL( SCIPsetIncludeBenders(scip->set, benders) );

   return SCIP_OKAY;
}

/** creates a Benders' decomposition and includes it in SCIP with all non-fundamental callbacks set to NULL
 *
 *  If needed, the non-fundamental callbacks can be added afterwards via setter functions SCIPsetBendersCopy(),
 *  SCIPsetBendersFree(), SCIPsetBendersInity(), SCIPsetBendersExit(), SCIPsetBendersInitsol(), SCIPsetBendersExitsol(),
 *  SCIPsetBendersFarkas().
 *
 *  To use the Benders' decomposition for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeBenders() instead
 */
SCIP_RETCODE SCIPincludeBendersBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS**        bendersptr,         /**< reference to a benders, or NULL */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   SCIP_BENDERS* benders;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBendersBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether Benders' decomposition is already present */
   if( SCIPfindBenders(scip, name) != NULL )
   {
      SCIPerrorMessage("benders <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbendersCreate(&benders, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         cutlp, cutpseudo, cutrelax, shareauxvars, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, bendersgetvar,
         benderscreatesub, NULL, NULL, NULL, NULL, NULL, bendersdata) );
   SCIP_CALL( SCIPsetIncludeBenders(scip->set, benders) );

   if( bendersptr != NULL )
      *bendersptr = benders;

   return SCIP_OKAY;
}

/** sets copy method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSCOPY((*benderscopy))     /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetCopy(benders, benderscopy);

   return SCIP_OKAY;
}

/** sets destructor method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREE((*bendersfree))     /**< destructor of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetFree(benders, bendersfree);

   return SCIP_OKAY;
}

/** sets initialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit))    /**< initialize Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetInit(benders, bendersinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit))    /**< deinitialize Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetExit(benders, bendersexit);

   return SCIP_OKAY;
}

/** sets presolving initialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersInitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))/**< presolving initialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersInitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetInitpre(benders, bendersinitpre);

   return SCIP_OKAY;
}

/** sets presolving deinitialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersExitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))/**< presolving deinitialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersExitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetExitpre(benders, bendersexitpre);

   return SCIP_OKAY;
}

/** sets solving process initialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetInitsol(benders, bendersinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetExitsol(benders, bendersexitsol);

   return SCIP_OKAY;
}

/** sets the method called prior to solving the subproblems for benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersPresubsolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve))/**< method called prior to solving the subproblems */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersPresubsolve", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetPresubsolve(benders, benderspresubsolve);

   return SCIP_OKAY;
}

/** sets the subproblem solving and freeing methods for Benders' decomposition
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersSolveAndFreesub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< solving method for a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub))/**< the subproblem freeing method for Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersSolveAndFreesub", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benderssolvesubconvex == NULL && benderssolvesub == NULL && bendersfreesub != NULL)
      || ((benderssolvesubconvex != NULL || benderssolvesub != NULL) && bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that if bendersFreesub%s is "
         "implemented at least one of bendersSolvesubconvex%s or bendersSolvesub%s are implemented, "
         "or if bendersFreesub%s is not implemented, then none are implented.\n", SCIPbendersGetName(benders),
         SCIPbendersGetName(benders), SCIPbendersGetName(benders), SCIPbendersGetName(benders),
         SCIPbendersGetName(benders));
      return SCIP_INVALIDCALL;
   }

   SCIPbendersSetSolvesubconvex(benders, benderssolvesubconvex);
   SCIPbendersSetSolvesub(benders, benderssolvesub);
   SCIPbendersSetFreesub(benders, bendersfreesub);

   return SCIP_OKAY;
}

/** sets the post solving methods for benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersPostsolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersPostsolve", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetPostsolve(benders, benderspostsolve);

   return SCIP_OKAY;
}

/** sets the subproblem comparison method for determining the solving order in Benders' decomposition
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersSubproblemComp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_SORTPTRCOMP((*benderssubcomp))  /**< a comparator for defining the solving order of the subproblems */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersSubproblemComp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetSubproblemComp(benders, benderssubcomp);

   return SCIP_OKAY;
}

/** returns the Benders' decomposition of the given name, or NULL if not existing */
SCIP_BENDERS* SCIPfindBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of Benders' decomposition */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindBenders(scip->set, name);
}

/** returns the array of currently available Benders' decomposition; active Benders' decomposition are in the first
 * slots of the array
 */
SCIP_BENDERS** SCIPgetBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortBenders(scip->set);

   return scip->set->benders;
}

/** returns the number of currently available Benders' decomposition */
int SCIPgetNBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nbenders;
}

/** returns the number of currently active Benders' decomposition */
int SCIPgetNActiveBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nactivebenders;
}

/** activates the Benders' decomposition to be used for the current problem
 *
 *  This method should be called during the problem creation stage for all pricers that are necessary to solve
 *  the problem model.
 *
 *  @note The Benders' decompositions are automatically deactivated when the problem is freed.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPactivateBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPactivateBenders", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersActivate(benders, scip->set, nsubproblems) );

   return SCIP_OKAY;
}

/** deactivates the Benders' decomposition
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPdeactivateBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdeactivateBenders", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersDeactivate(benders, scip->set) );

   return SCIP_OKAY;
}

/** sets the priority of a Benders' decomposition */
void SCIPsetBendersPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);

   SCIPbendersSetPriority(benders, scip->set, priority);
}

/** calls the exec method of Benders' decomposition to solve the subproblems
 *
 *  The checkint flag indicates whether integer feasibility can be assumed. If it is not assumed, i.e. checkint ==
 *  FALSE, then only the convex relaxations of the subproblems are solved. If integer feasibility is assumed, i.e.
 *  checkint == TRUE, then the convex relaxations and the full CIP are solved to generate Benders' cuts and check
 *  solution feasibility.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPsolveBendersSubproblems(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            auxviol,            /**< set to TRUE only if the solution is feasible but the aux vars are violated */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should the integer solution be checked by the subproblems */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveBendersSubproblems", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersExec(benders, scip->set, sol, result, infeasible, auxviol, type, checkint) );

   return SCIP_OKAY;
}

/** returns the master problem variable for the given subproblem variable
 *
 *  This function is used as part of the cut generation process.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetBendersMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_VAR*             var,                /**< the subproblem variable */
   SCIP_VAR**            mappedvar           /**< pointer to store the master variable that var is mapped to */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetBendersMasterVar", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersGetVar(benders, scip->set, var, mappedvar, -1) );

   return SCIP_OKAY;
}

/** returns the subproblem problem variable for the given master variable
 *
 *  This function is used as part of the cut generation process.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetBendersSubproblemVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_VAR*             var,                /**< the master variable */
   SCIP_VAR**            mappedvar,          /**< pointer to store the subproblem variable that var is mapped to */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetBendersSubproblemVar", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersGetVar(benders, scip->set, var, mappedvar, probnumber) );

   return SCIP_OKAY;
}

/** returns the number of subproblems that are stored in the given Benders' decomposition
 *
 *  @return the number of subproblems in the Benders' decomposition
 */
int SCIPgetBendersNSubproblems(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);

   return SCIPbendersGetNSubproblems(benders);
}

/** registers the Benders' decomposition subproblem with the Benders' decomposition struct.
 *
 *  If a custom subproblem solving method is used and no internal cut generation methods will be employed, then the
 *  subproblem parameter can be set to NULL. By setting subproblem to NULL will inform the Benders' decomposition core
 *  that a custom solving method is used. This will ensure that no internal solving methods are invoked during the
 *  solution process.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 */
SCIP_RETCODE SCIPaddBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP*                 subproblem          /**< Benders' decomposition subproblem, can be NULL */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddBendersSubproblem", FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersAddSubproblem(benders, subproblem) );

   return SCIP_OKAY;
}

/** calls the generic subproblem setup method for a Benders' decomposition subproblem
 *
 *  This is called if the user requires to solve the Benders' decomposition subproblem separately from the main Benders'
 *  solving loop. This could be in the case of enhancement techniques.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPsetupBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   SCIP_SOL*             sol,                /**< primal solution used to setup the problem, NULL for LP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetupBendersSubproblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersSetupSubproblem(benders, scip->set, sol, probnumber, type) );

   return SCIP_OKAY;
}

/** calls the solving method for a single Benders' decomposition subproblem
 *
 *  The method either calls the users solve subproblem method or calls the generic method. In the case of the generic
 *  method, the user must set up the subproblem prior to calling this method.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPsolveBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL for the current LP/Pseudo solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_Bool             solvecip,           /**< directly solve the CIP subproblem */
   SCIP_Real*            objective           /**< the objective function value of the subproblem, can be NULL */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveBendersSubproblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

  SCIP_CALL( SCIPbendersSolveSubproblem(benders, scip->set, sol, probnumber, infeasible, solvecip, objective) );

   return SCIP_OKAY;
}

/** frees the subproblem after calling the solve subproblem method
 *
 *  This will either call the user defined free
 *  subproblem callback for Benders' decomposition or the default freeing methods. In the default case, if the
 *  subproblem is an LP, then SCIPendProbing is called. If the subproblem is a MIP, then SCIPfreeTransform is called.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPfreeBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeBendersSubproblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersFreeSubproblem(benders, scip->set, probnumber) );

   return SCIP_OKAY;
}

/** checks the optimality of a Benders' decomposition subproblem by comparing the objective function value against the
 *  value of the corresponding auxiliary variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called in any SCIP stage; however, it is necessary that the Benders' auxiliary variables
 *  exist. If they don't exist, then the optimal flag is returned as FALSE.
 *
 *  @pre This method can be called if requested subproblem is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcheckBendersSubproblemOptimality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL for the current LP solution */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_Bool*            optimal             /**< flag to indicate whether the current subproblem is optimal for the master */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   (*optimal) = FALSE;

   if( SCIPbendersGetAuxiliaryVar(benders, probnumber) == NULL )
   {
      SCIPinfoMessage(scip, NULL, "Benders' decomposition: The auxiliary variable for subproblem <%d> doesn't exist. "
         "SCIPcheckBendersSubproblemOptimality can not be currently called at stage <%d>.\n", probnumber,
         SCIPgetStage(scip));
      SCIPinfoMessage(scip, NULL, "  The optimal flag will be returned as FALSE.\n");

      return SCIP_OKAY;
   }

   if( SCIPbendersSubproblem(benders, probnumber) != NULL )
   {
      SCIP_CALL( SCIPcheckStage(SCIPbendersSubproblem(benders, probnumber), "SCIPcheckBendersSubproblemOptimality",
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   }

   (*optimal) = SCIPbendersSubproblemIsOptimal(benders, scip->set, sol, probnumber);

   return SCIP_OKAY;
}

/** returns the value of the auxiliary variable for a given subproblem
 *
 *  @return the value of the auxiliary variable for the given subproblem
 */
SCIP_Real SCIPgetBendersAuxiliaryVarVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL for the current LP solution */
   int                   probnumber          /**< the number of the pricing problem */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return SCIPbendersGetAuxiliaryVarVal(benders, scip->set, sol, probnumber);
}

/** solves an independent subproblem to identify its lower bound and updates the lower bound of the corresponding
 *  auxiliary variable
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcomputeBendersSubproblemLowerbound(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem to be evaluated */
   SCIP_Real*            lowerbound,         /**< the lowerbound for the subproblem */
   SCIP_Bool*            infeasible          /**< was the subproblem found to be infeasible? */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcomputeBendersSubproblemLowerbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersComputeSubproblemLowerbound(benders, scip->set, probnumber, lowerbound, infeasible) );

   return SCIP_OKAY;
}

/** Merges a subproblem into the master problem. This process just adds a copy of the subproblem variables and
 *  constraints to the master problem, but keeps the subproblem stored in the Benders' decomposition data structure.
 *  The reason for keeping the subproblem available is for when it is queried for solutions after the problem is solved.
 *
 *  Once the subproblem is merged into the master problem, then the subproblem is flagged as disabled. This means that
 *  it will not be solved in the subsequent subproblem solving loops.
 *
 *  The associated auxiliary variables are kept in the master problem. The objective function of the merged subproblem
 *  is added as an underestimator constraint.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPmergeBendersSubproblemIntoMaster(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of subproblem variables corresponding
                                              *   to the newly created master variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of subproblem constraints to the
                                              *   corresponding newly created constraints, or NULL */
   int                   probnumber          /**< the number of the subproblem that will be merged into the master problem*/
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPmergeBendersSubproblemIntoMaster", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersMergeSubproblemIntoMaster(benders, scip->set, varmap, consmap, probnumber) );

   return SCIP_OKAY;
}

/** applies a Benders' decomposition to the selected decomposition from the decomposition store
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPapplyBendersDecomposition(
   SCIP*                 scip,               /**< the SCIP data structure */
   int                   decompindex         /**< the index of the decomposition that will be applied */
   )
{
   SCIP_BENDERS* benders;

   assert(scip != NULL);
   assert(scip->decompstore != NULL);
   assert(SCIPdecompstoreGetNOrigDecomps(scip->decompstore) > 0);
   assert(decompindex >= 0 && decompindex < SCIPdecompstoreGetNOrigDecomps(scip->decompstore));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPapplyBendersDecomposition", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* if there already exists an active Benders' decomposition, then default decomposition is not applied. */
   if( SCIPgetNActiveBenders(scip) > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "A Benders' decomposition already exists. The default Benders' decomposition will not be applied to the stored decomposition.\n");
      return SCIP_OKAY;
   }

   /* retrieving the default Benders' decomposition plugin */
   benders = SCIPfindBenders(scip, "default");

   /* if the default Benders' decomposition plugin doesn't exist, then this will result in an error */
   if( benders == NULL )
   {
      SCIPerrorMessage("The default Benders' decomposition plugin is required to apply Benders' decomposition using the input decomposition.");
      return SCIP_ERROR;
   }

   /* applying the Benders' decomposition. If SCIP is in the PROBLEM stage, then the auxiliary variables don't need to
    * be added. However, in any other stage, then the auxiliary variables must be added to the problem.
    */
   SCIP_CALL( SCIPbendersApplyDecomposition(benders, scip->set,
         SCIPdecompstoreGetOrigDecomps(scip->decompstore)[decompindex]) );

   return SCIP_OKAY;
}

/** creates a Benders' cut algorithms and includes it in the associated Benders' decomposition
 *
 *  This should be called from the SCIPincludeBendersXyz for the associated Benders' decomposition. It is only possible
 *  to include a Benders' cut algorithm if a Benders' decomposition has already been included
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note method has all Benders' decomposition callbacks as arguments and is thus changed every time a new callback is
 *        added in future releases; consider using SCIPincludeBendersBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeBenderscut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name,               /**< name of Benders' decomposition cuts */
   const char*           desc,               /**< description of Benders' decomposition cuts */
   int                   priority,           /**< priority of the Benders' decomposition cuts */
   SCIP_Bool             islpcut,            /**< indicates whether the cut is generated from the LP solution */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy)),/**< copy method of Benders' decomposition cuts or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree)),/**< destructor of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit)),/**< initialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit)),/**< deinitialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol)),/**< solving process initialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol)),/**< solving process deinitialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< execution method of Benders' decomposition cuts */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' decomposition cuts data */
   )
{
   SCIP_BENDERSCUT* benderscut;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBenderscut", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether pricer is already present */
   if( SCIPfindBenderscut(benders, name) != NULL )
   {
      SCIPerrorMessage("benders <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbenderscutCreate(benders, &benderscut, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         islpcut, benderscutcopy, benderscutfree, benderscutinit, benderscutexit,
         benderscutinitsol, benderscutexitsol, benderscutexec, benderscutdata) );
   SCIP_CALL( SCIPbendersIncludeBenderscut(benders, scip->set, benderscut) );

   return SCIP_OKAY;
}

/** creates a Benders' cut and includes it an associated Benders' decomposition with all non-fundamental callbacks set to NULL
 *
 *  If needed, the non-fundamental callbacks can be added afterwards via setter functions SCIPsetBenderscutCopy(),
 *  SCIPsetBenderscutFree(), SCIPsetBenderscutInit(), SCIPsetBenderscutExit(), SCIPsetBenderscutInitsol(),
 *  SCIPsetBenderscutExitsol().
 *
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeBenders() instead
 */
SCIP_RETCODE SCIPincludeBenderscutBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT**     benderscutptr,      /**< reference to a Benders' decomposition cut, or NULL */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             islpcut,            /**< indicates whether the cut is generated from the LP solution */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< the execution method of the Benders' cut algorithm */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' cut data */
   )
{
   SCIP_BENDERSCUT* benderscut;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBenderscutBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether Benders' decomposition cut is already present */
   if( SCIPfindBenderscut(benders, name) != NULL )
   {
      SCIPerrorMessage("Benders' cut <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbenderscutCreate(benders, &benderscut, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc,
         priority, islpcut, NULL, NULL, NULL, NULL, NULL, NULL, benderscutexec, benderscutdata) );
   SCIP_CALL( SCIPbendersIncludeBenderscut(benders, scip->set, benderscut) );

   if( benderscutptr != NULL )
      *benderscutptr = benderscut;

   return SCIP_OKAY;
}

/** sets copy method of Benders' decomposition cut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy))/**< copy method of benderscut or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetCopy(benderscut, benderscutcopy);

   return SCIP_OKAY;
}

/** sets destructor method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree))/**< destructor of benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetFree(benderscut, benderscutfree);

   return SCIP_OKAY;
}

/** sets initialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit))/**< initialize benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetInit(benderscut, benderscutinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit))/**< deinitialize benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetExit(benderscut, benderscutexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol))/**< solving process initialization method of benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetInitsol(benderscut, benderscutinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol))/**< solving process deinitialization method of benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetExitsol(benderscut, benderscutexitsol);

   return SCIP_OKAY;
}

/** sets the priority of a Benders' decomposition cut algorithm
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   SCIP_BENDERS** benders;
   int nbenders;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benderscut != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutPriority", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPbenderscutSetPriority(benderscut, priority);

   /* DIRTY: This is not a good fix */
   /* Changing the priority of one Benders' cut in a Benders' decomposition requires all Benders' cuts to be set to
    * unsorted. This is a fix that is not very nice, but it does the job */
   benders = SCIPgetBenders(scip);
   nbenders = SCIPgetNBenders(scip);
   for( i = 0; i < nbenders; i++ )
      SCIPbendersSetBenderscutsSorted(benders[i], FALSE);

   return SCIP_OKAY;
}

/** adds the generated cuts to the Benders' cut storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPstoreBendersCut(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_VAR**            vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real*            vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real             lhs,                /**< the left hand side of the cut */
   SCIP_Real             rhs,                /**< the right hand side of the cut */
   int                   nvars               /**< the number of variables with non-zero coefficients in the cut */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstoreBendersCut", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersStoreCut(benders, scip->set, vars, vals, lhs, rhs, nvars) );

   return SCIP_OKAY;
}

/** creates a constraint in the input SCIP instance that corresponds to the given vars and vals arrays */
static
SCIP_RETCODE createAndApplyStoredBendersCut(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_VAR**            vars,               /**< the variables from the source constraint */
   SCIP_Real*            vals,               /**< the coefficients of the variables in the source constriant */
   SCIP_Real             lhs,                /**< the LHS of the source constraint */
   SCIP_Real             rhs,                /**< the RHS of the source constraint */
   int                   nvars,              /**< the number of variables in the source constraint */
   int                   consindex           /**< the store index for the constraint */
   )
{
   SCIP_CONS* cons;
   char consname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   /* setting the name of the transferred cut */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "transferredbenderscut_%d", consindex);

   /* creating a linear constraint given the input parameters */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, nvars, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPsetConsRemovable(scip, cons, TRUE) );

   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** applies the Benders' decomposition cuts in storage to the input SCIP instance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPapplyBendersStoredCuts(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int naddedcuts;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   /* retrieving the number of stored Benders' cuts */
   naddedcuts =  SCIPbendersGetNStoredCuts(benders);

   /* looping over all added cuts to construct the cut for the input SCIP instance */
   for( i = 0; i < naddedcuts; i++ )
   {
      /* collecting the variable information from the constraint */
      SCIP_CALL( SCIPbendersGetStoredCutData(benders, i, &vars, &vals, &lhs, &rhs, &nvars) );

      if( nvars > 0 )
      {
         /* create and apply the cut to be transferred from the sub SCIP to the source SCIP */
         SCIP_CALL( createAndApplyStoredBendersCut(scip, vars, vals, lhs, rhs, nvars, i) );
      }
   }

   return SCIP_OKAY;
}
