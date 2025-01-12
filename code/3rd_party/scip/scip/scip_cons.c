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

/**@file   scip_cons.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for constraint handler plugins and constraints
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

#include "scip/cons.h"
#include "scip/debug.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/set.h"
#include "scip/struct_cons.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPmarkConsPropagate

/** creates a constraint handler and includes it in SCIP.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note method has all constraint handler callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludeConshdlrBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of constraint handler */
   const char*           desc,               /**< description of constraint handler */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   int                   enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int                   chckpriority,       /**< priority of the constraint handler for checking feasibility (and propagation) */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int                   eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   SCIP_Bool             delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   SCIP_PROPTIMING       proptiming,         /**< positions in the node solving loop where propagation method of constraint handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the constraint handler's presolving method */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy)),  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   SCIP_DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   SCIP_DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   SCIP_DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   SCIP_DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   SCIP_DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   SCIP_DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   SCIP_DECL_CONSENFORELAX ((*consenforelax)), /**< enforcing constraints for relaxation solutions */
   SCIP_DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   SCIP_DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   SCIP_DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   SCIP_DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   SCIP_DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   SCIP_DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   SCIP_DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   SCIP_DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   SCIP_DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   SCIP_DECL_CONSDELVARS ((*consdelvars)),   /**< variable deletion method */
   SCIP_DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   SCIP_DECL_CONSCOPY    ((*conscopy)),      /**< constraint copying method */
   SCIP_DECL_CONSPARSE   ((*consparse)),     /**< constraint parsing method */
   SCIP_DECL_CONSGETVARS ((*consgetvars)),   /**< constraint get variables method */
   SCIP_DECL_CONSGETNVARS((*consgetnvars)),  /**< constraint get number of variable method */
   SCIP_DECL_CONSGETDIVEBDCHGS((*consgetdivebdchgs)), /**< constraint handler diving solution enforcement method */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeConshdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether constraint handler is already present */
   if( SCIPfindConshdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("constraint handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconshdlrCreate(&conshdlr, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, sepapriority, enfopriority, chckpriority, sepafreq, propfreq, eagerfreq, maxprerounds,
         delaysepa, delayprop, needscons, proptiming, presoltiming, conshdlrcopy,
         consfree, consinit, consexit, consinitpre, consexitpre, consinitsol, consexitsol,
         consdelete, constrans, consinitlp, conssepalp, conssepasol, consenfolp, consenforelax, consenfops, conscheck, consprop,
         conspresol, consresprop, conslock, consactive, consdeactive, consenable, consdisable, consdelvars, consprint,
         conscopy, consparse, consgetvars, consgetnvars, consgetdivebdchgs, conshdlrdata) );
   SCIP_CALL( SCIPsetIncludeConshdlr(scip->set, conshdlr) );

   return SCIP_OKAY;
}

/** creates a constraint handler and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetConshdlrInit(), SCIPsetConshdlrExit(),
 *  SCIPsetConshdlrCopy(), SCIPsetConshdlrFree(), SCIPsetConshdlrInitsol(), SCIPsetConshdlrExitsol(),
 *  SCIPsetConshdlrInitpre(), SCIPsetConshdlrExitpre(), SCIPsetConshdlrPresol(), SCIPsetConshdlrDelete(),
 *  SCIPsetConshdlrDelvars(), SCIPsetConshdlrInitlp(), SCIPsetConshdlrActive(), SCIPsetConshdlrDeactive(),
 *  SCIPsetConshdlrEnable(), SCIPsetConshdlrDisable(), SCIPsetConshdlrResprop(), SCIPsetConshdlrTrans(),
 *  SCIPsetConshdlrPrint(), SCIPsetConshdlrParse(), SCIPsetConshdlrGetVars(), SCIPsetConshdlrGetNVars(), and
 *  SCIPsetConshdlrGetDiveBdChgs().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeConshdlr() instead
 */
SCIP_RETCODE SCIPincludeConshdlrBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR**       conshdlrptr,        /**< reference to a constraint handler pointer, or NULL */
   const char*           name,               /**< name of constraint handler */
   const char*           desc,               /**< description of constraint handler */
   int                   enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int                   chckpriority,       /**< priority of the constraint handler for checking feasibility (and propagation) */
   int                   eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   SCIP_Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   SCIP_DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   SCIP_DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   SCIP_DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   SCIP_DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeConshdlrBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether constraint handler is already present */
   if( SCIPfindConshdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("constraint handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconshdlrCreate(&conshdlr, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, 0, enfopriority, chckpriority, -1, -1, eagerfreq, 0,
         FALSE, FALSE, needscons,
         SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_ALWAYS,
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL, consenfolp, NULL, consenfops, conscheck, NULL,
         NULL, NULL, conslock, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL, conshdlrdata) );
   SCIP_CALL( SCIPsetIncludeConshdlr(scip->set, conshdlr) );

   if( conshdlrptr != NULL )
      *conshdlrptr = conshdlr;

   return SCIP_OKAY;
}

/** sets all separation related callbacks/parameters of the constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   SCIP_Bool             delaysepa           /**< should separation method be delayed, if other separators found cuts? */
   )
{
   int oldsepapriority;
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrSepa", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   oldsepapriority = SCIPconshdlrGetSepaPriority(conshdlr);
   SCIPconshdlrSetSepa(conshdlr, conssepalp, conssepasol, sepafreq, sepapriority, delaysepa);

   /* change the position of the constraint handler in the constraint handler array w.r.t. its new sepa priority */
   if( oldsepapriority != sepapriority )
      SCIPsetReinsertConshdlrSepaPrio(scip->set, conshdlr, oldsepapriority);

   name = SCIPconshdlrGetName(conshdlr);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/sepafreq", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, sepafreq) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/delaysepa", name);
   SCIP_CALL( SCIPsetSetDefaultBoolParam(scip->set, paramname, delaysepa) );

   return SCIP_OKAY;
}

/** sets both the propagation callback and the propagation frequency of the constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       proptiming          /**< positions in the node solving loop where propagation should be executed */
   )
{
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrProp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetProp(conshdlr, consprop, propfreq, delayprop, proptiming);

   name = SCIPconshdlrGetName(conshdlr);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/propfreq", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, propfreq) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/proptiming", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, (int) proptiming) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/delayprop", name);
   SCIP_CALL( SCIPsetSetDefaultBoolParam(scip->set, paramname, delayprop) );

   return SCIP_OKAY;
}

/** sets relaxation enforcement method of the constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrEnforelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSENFORELAX ((*consenforelax)) /**< enforcement method for relaxation solution of constraint handler (might be NULL) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrEnforelax", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetEnforelax(conshdlr, consenforelax);

   return SCIP_OKAY;
}

/** sets copy method of both the constraint handler and each associated constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy)),  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSCOPY    ((*conscopy))       /**< constraint copying method */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetCopy(conshdlr, conshdlrcopy, conscopy);

   return SCIP_OKAY;
}

/** sets destructor method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSFREE    ((*consfree))       /**< destructor of constraint handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetFree(conshdlr, consfree);

   return SCIP_OKAY;
}

/** sets initialization method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit))       /**< initialize constraint handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetInit(conshdlr, consinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit))       /**< deinitialize constraint handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetExit(conshdlr, consexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINITSOL((*consinitsol))     /**< solving process initialization method of constraint handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetInitsol(conshdlr, consinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol))/**< solving process deinitialization method of constraint handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetExitsol(conshdlr, consexitsol);

   return SCIP_OKAY;
}

/** sets preprocessing initialization method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrInitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINITPRE((*consinitpre))     /**< preprocessing initialization method of constraint handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrInitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetInitpre(conshdlr, consinitpre);

   return SCIP_OKAY;
}

/** sets preprocessing deinitialization method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrExitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSEXITPRE((*consexitpre))     /**< preprocessing deinitialization method of constraint handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrExitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(conshdlr != NULL);

   SCIPconshdlrSetExitpre(conshdlr, consexitpre);

   return SCIP_OKAY;
}

/** sets presolving method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method of constraint handler */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the constraint handler's presolving method */
   )
{
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrPresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconshdlrSetPresol(conshdlr, conspresol, maxprerounds, presoltiming) );

   name = SCIPconshdlrGetName(conshdlr);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, maxprerounds) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/%s/presoltiming", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, (int) presoltiming) );

   return SCIP_OKAY;
}

/** sets method of constraint handler to free specific constraint data
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete))     /**< free specific constraint data */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrDelete", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetDelete(conshdlr, consdelete);

   return SCIP_OKAY;
}

/** sets method of constraint handler to transform constraint data into data belonging to the transformed problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrTrans(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSTRANS   ((*constrans))      /**< transform constraint data into data belonging to the transformed problem */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrTrans", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetTrans(conshdlr, constrans);

   return SCIP_OKAY;
}

/** sets method of constraint handler to initialize LP with relaxations of "initial" constraints
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrInitlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSINITLP  ((*consinitlp))     /**< initialize LP with relaxations of "initial" constraints */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrInitlp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetInitlp(conshdlr, consinitlp);

   return SCIP_OKAY;
}

/** sets propagation conflict resolving method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrResprop(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSRESPROP ((*consresprop))    /**< propagation conflict resolving method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrResprop", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetResprop(conshdlr, consresprop);

   return SCIP_OKAY;
}

/** sets activation notification method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrActive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSACTIVE  ((*consactive))     /**< activation notification method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrActive", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetActive(conshdlr, consactive);

   return SCIP_OKAY;
}

/** sets deactivation notification method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrDeactive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDEACTIVE((*consdeactive))   /**< deactivation notification method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrDeactive", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetDeactive(conshdlr, consdeactive);

   return SCIP_OKAY;
}

/** sets enabling notification method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrEnable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSENABLE  ((*consenable))     /**< enabling notification method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrEnable", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetEnable(conshdlr, consenable);

   return SCIP_OKAY;
}

/** sets disabling notification method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrDisable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDISABLE ((*consdisable))    /**< disabling notification method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrDisable", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetDisable(conshdlr, consdisable);

   return SCIP_OKAY;
}

/** sets variable deletion method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrDelvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSDELVARS ((*consdelvars))    /**< variable deletion method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrDelvars", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetDelvars(conshdlr, consdelvars);

   return SCIP_OKAY;
}

/** sets constraint display method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPRINT   ((*consprint))      /**< constraint display method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrPrint", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetPrint(conshdlr, consprint);

   return SCIP_OKAY;
}

/** sets constraint parsing method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrParse(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPARSE   ((*consparse))      /**< constraint parsing method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrParse", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetParse(conshdlr, consparse);

   return SCIP_OKAY;
}

/** sets constraint variable getter method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrGetVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSGETVARS ((*consgetvars))    /**< constraint variable getter method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrGetVars", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetGetVars(conshdlr, consgetvars);

   return SCIP_OKAY;
}

/** sets constraint variable number getter method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrGetNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSGETNVARS((*consgetnvars))   /**< constraint variable number getter method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrGetNVars", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetGetNVars(conshdlr, consgetnvars);

   return SCIP_OKAY;
}

/** sets diving bound change method of constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetConshdlrGetDiveBdChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSGETDIVEBDCHGS((*consgetdivebdchgs)) /**< constraint handler diving solution enforcement method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConshdlrGetDiveBdChgs", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconshdlrSetGetDiveBdChgs(conshdlr, consgetdivebdchgs);

   return SCIP_OKAY;
}
/** returns the constraint handler of the given name, or NULL if not existing */
/** returns the constraint handler of the given name, or NULL if not existing */
SCIP_CONSHDLR* SCIPfindConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint handler */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindConshdlr(scip->set, name);
}

/** returns the array of currently available constraint handlers */
SCIP_CONSHDLR** SCIPgetConshdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->conshdlrs;
}

/** returns the number of currently available constraint handlers */
int SCIPgetNConshdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nconshdlrs;
}

/** creates and captures a constraint of the given constraint handler
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP_CONSDATA*        consdata,           /**< data for this specific constraint */
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
   assert(cons != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPconsCreate(cons, scip->mem->probmem, scip->set, name, conshdlr, consdata,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, TRUE, TRUE) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_EXITSOLVE:
      SCIP_CALL( SCIPconsCreate(cons, scip->mem->probmem, scip->set, name, conshdlr, consdata,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, FALSE, TRUE) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** parses constraint information (in cip format) out of a string; if the parsing process was successful a constraint is
 *  creates and captures;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 */
SCIP_RETCODE SCIPparseCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to store constraint */
   const char*           str,                /**< string to parse for constraint */
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
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool*            success             /**< pointer to store if the paring process was successful */
   )
{
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPparseCons", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsParse(cons, scip->set, scip->messagehdlr, str,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, success) );

   return SCIP_OKAY;
}

/** increases usage counter of constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcaptureCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to capture */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcaptureCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( cons->scip == scip );

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** decreases usage counter of constraint, if the usage pointer reaches zero the constraint gets freed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note the pointer of the constraint will be NULLed
 */
SCIP_RETCODE SCIPreleaseCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons                /**< pointer to constraint */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPreleaseCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPconsRelease(cons, scip->mem->probmem, scip->set) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      if( SCIPconsIsOriginal(*cons) && (*cons)->nuses == 1 && (*cons)->transorigcons != NULL )
      {
         SCIPerrorMessage("cannot release last use of original constraint while an associated transformed constraint exists\n");
         return SCIP_INVALIDCALL;
      }
      SCIP_CALL( SCIPconsRelease(cons, scip->mem->probmem, scip->set) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** change constraint name
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note to get the current name of a constraint, use SCIPconsGetName() from pub_cons.h
 */
SCIP_RETCODE SCIPchgConsName(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   const char*           name                /**< new name of constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgConsName", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE , FALSE, FALSE, FALSE) );

   assert( cons->scip == scip );

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("constraint names can only be changed in problem creation stage\n");
      SCIPABORT();
      return SCIP_INVALIDCALL; /*lint !e527*/
   }

   /* remove constraint's name from the namespace if the constraint was already added */
   if( SCIPconsIsAdded(cons) )
   {
      SCIP_CALL( SCIPprobRemoveConsName(scip->origprob, cons) );
   }

   /* change constraint name */
   SCIP_CALL( SCIPconsChgName(cons, SCIPblkmem(scip), name) );

   /* add constraint's name to the namespace if the constraint was already added */
   if( SCIPconsIsAdded(cons) )
   {
      SCIP_CALL( SCIPprobAddConsName(scip->origprob, cons) );
   }

   return SCIP_OKAY;
}

/** sets the initial flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsInitial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             initial             /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsInitial", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetInitial(cons, scip->set, scip->stat, initial) );

   return SCIP_OKAY;
}

/** sets the separate flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsSeparated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             separate            /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsSeparated", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetSeparated(cons, scip->set, separate) );

   return SCIP_OKAY;
}

/** sets the enforce flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsEnforced(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             enforce             /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsEnforced", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetEnforced(cons, scip->set, enforce) );

   return SCIP_OKAY;
}

/** sets the check flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsChecked(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             check               /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsChecked", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetChecked(cons, scip->set, check) );

   return SCIP_OKAY;
}

/** sets the propagate flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsPropagated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             propagate           /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsPropagated", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetPropagated(cons, scip->set, propagate) );

   return SCIP_OKAY;
}

/** sets the local flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             local               /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsLocal", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconsSetLocal(cons, local);

   return SCIP_OKAY;
}

/** sets the modifiable flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPsetConsModifiable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             modifiable          /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsModifiable", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   SCIPconsSetModifiable(cons, modifiable);

   return SCIP_OKAY;
}

/** sets the dynamic flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsDynamic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             dynamic             /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsDynamic", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconsSetDynamic(cons, dynamic);

   return SCIP_OKAY;
}

/** sets the removable flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsRemovable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             removable           /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsRemovable", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconsSetRemovable(cons, removable);

   return SCIP_OKAY;
}

/** sets the stickingatnode flag of the given constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetConsStickingAtNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             stickingatnode      /**< new value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConsStickingAtNode", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconsSetStickingAtNode(cons, stickingatnode);

   return SCIP_OKAY;
}

/** updates the flags of the first constraint according to the ones of the second constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPupdateConsFlags(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1               /**< constraint that should be deleted */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateConsFlags", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPconsIsInitial(cons1) )
   {
      SCIP_CALL( SCIPsetConsInitial(scip, cons0, TRUE) );
   }
   if( SCIPconsIsSeparated(cons1) )
   {
      SCIP_CALL( SCIPsetConsSeparated(scip, cons0, TRUE) );
   }
   if( SCIPconsIsEnforced(cons1) )
   {
      SCIP_CALL( SCIPsetConsEnforced(scip, cons0, TRUE) );
   }
   if( SCIPconsIsChecked(cons1) )
   {
      SCIP_CALL( SCIPsetConsChecked(scip, cons0, TRUE) );
   }
   if( SCIPconsIsPropagated(cons1) )
   {
      SCIP_CALL( SCIPsetConsPropagated(scip, cons0, TRUE) );
   }
   if( !SCIPconsIsDynamic(cons1) )
   {
      SCIP_CALL( SCIPsetConsDynamic(scip, cons0, FALSE) );
   }
   if( !SCIPconsIsRemovable(cons1) )
   {
      SCIP_CALL( SCIPsetConsRemovable(scip, cons0, FALSE) );
   }
   if( SCIPconsIsStickingAtNode(cons1) )
   {
      SCIP_CALL( SCIPsetConsStickingAtNode(scip, cons0, TRUE) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPtransformCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get/create transformed constraint for */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(transcons != NULL);
   assert(cons->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtransformCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPconsIsTransformed(cons) )
   {
      *transcons = cons;
      SCIPconsCapture(*transcons);
   }
   else
   {
      SCIP_CALL( SCIPconsTransform(cons, scip->mem->probmem, scip->set, transcons) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed constraints for an array of constraints;
 *  if a constraint in the array is not yet transformed, a new transformed constraint for this constraint is created;
 *  it is possible to call this method with conss == transconss
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPtransformConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nconss,             /**< number of constraints to get/create transformed constraints for */
   SCIP_CONS**           conss,              /**< array with constraints to get/create transformed constraints for */
   SCIP_CONS**           transconss          /**< array to store the transformed constraints */
   )
{
   int c;

   assert(nconss == 0 || conss != NULL);
   assert(nconss == 0 || transconss != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtransformConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsTransformed(conss[c]) )
      {
         transconss[c] = conss[c];
         SCIPconsCapture(transconss[c]);
      }
      else
      {
         SCIP_CALL( SCIPconsTransform(conss[c], scip->mem->probmem, scip->set, &transconss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed constraint of a given constraint;
 *  returns NULL as transcons, if transformed constraint is not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
SCIP_RETCODE SCIPgetTransformedCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get the transformed constraint for */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(transcons != NULL);
   assert(cons->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetTransformedCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPconsIsTransformed(cons) )
      *transcons = cons;
   else
      *transcons = SCIPconsGetTransformed(cons);

   return SCIP_OKAY;
}

/** gets corresponding transformed constraints for an array of constraints;
 *  stores NULL in a transconss slot, if the transformed constraint is not yet existing;
 *  it is possible to call this method with conss == transconss, but remember that constraints that are not
 *  yet transformed will be replaced with NULL
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
SCIP_RETCODE SCIPgetTransformedConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nconss,             /**< number of constraints to get the transformed constraints for */
   SCIP_CONS**           conss,              /**< constraints to get the transformed constraints for */
   SCIP_CONS**           transconss          /**< array to store the transformed constraints */
   )
{
   int c;

   assert(nconss == 0 || conss != NULL);
   assert(nconss == 0 || transconss != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetTransformedConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsTransformed(conss[c]) )
         transconss[c] = conss[c];
      else
         transconss[c] = SCIPconsGetTransformed(conss[c]);
   }

   return SCIP_OKAY;
}

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPaddConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             deltaage            /**< value to add to the constraint's age */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConsAge", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsAddAge(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob, deltaage, scip->reopt) );

   return SCIP_OKAY;
}

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPincConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPincConsAge", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsIncAge(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->reopt) );

   return SCIP_OKAY;
}

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPresetConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPresetConsAge", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsResetAge(cons, scip->set) );

   return SCIP_OKAY;
}

/** enables constraint's separation, propagation, and enforcing capabilities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPenableCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPenableCons", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnable(cons, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** disables constraint's separation, propagation, and enforcing capabilities, s.t. the constraint is not propagated,
 *  separated, and enforced anymore until it is enabled again with a call to SCIPenableCons();
 *  in contrast to SCIPdelConsLocal() and SCIPdelConsNode(), the disabling is not associated to a node in the tree and
 *  does not consume memory; therefore, the constraint is neither automatically enabled on leaving the node nor
 *  automatically disabled again on entering the node again;
 *  note that the constraints enforcing capabilities are necessary for the solution's feasibility, if the constraint
 *  is a model constraint; that means, you must be sure that the constraint cannot be violated in the current subtree,
 *  and you have to enable it again manually by calling SCIPenableCons(), if this subtree is left (e.g. by using
 *  an appropriate event handler that watches the corresponding variables' domain changes)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPdisableCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdisableCons", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsDisable(cons, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** enables constraint's separation capabilities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPenableConsSeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPenableConsSeparation", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnableSeparation(cons, scip->set) );

   return SCIP_OKAY;
}

/** disables constraint's separation capabilities s.t. the constraint is not separated anymore until the separation
 *  is enabled again with a call to SCIPenableConsSeparation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPdisableConsSeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdisableConsSeparation", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsDisableSeparation(cons, scip->set) );

   return SCIP_OKAY;
}

/** enables constraint's propagation capabilities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPenableConsPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPenableConsPropagation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnablePropagation(cons, scip->set) );

   return SCIP_OKAY;
}

/** disables constraint's propagation capabilities s.t. the constraint is not propagated anymore until the propagation
 *  is enabled again with a call to SCIPenableConsPropagation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPdisableConsPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdisableConsPropagation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsDisablePropagation(cons, scip->set) );

   return SCIP_OKAY;
}

#undef SCIPmarkConsPropagate

/** marks constraint to be propagated
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note if a constraint is marked to be propagated, the age of the constraint will be ignored for propagation
 */
SCIP_RETCODE SCIPmarkConsPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPmarkConsPropagate", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsMarkPropagate(cons, scip->set) );

   assert(!SCIPconsIsEnabled(cons) || SCIPconsIsMarkedPropagate(cons));

   return SCIP_OKAY;
}

/** unmarks the constraint to be propagated
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPunmarkConsPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPunmarkConsPropagate", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsUnmarkPropagate(cons, scip->set) );

   assert(!SCIPconsIsEnabled(cons) || !SCIPconsIsMarkedPropagate(cons));

   return SCIP_OKAY;
}

/** adds given values to lock status of type @p locktype of the constraint and updates the rounding locks of the involved variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPaddConsLocksType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_LOCKTYPE         locktype,           /**< type of the variable locks */
   int                   nlockspos,          /**< increase in number of rounding locks for constraint */
   int                   nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConsLocksType", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPconsAddLocks(cons, scip->set, locktype, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note This methods always adds locks of type model
 */
SCIP_RETCODE SCIPaddConsLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nlockspos,          /**< increase in number of rounding locks for constraint */
   int                   nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConsLocks", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPaddConsLocksType(scip, cons, SCIP_LOCKTYPE_MODEL, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}

/** checks single constraint for feasibility of the given solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcheckCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcheckCons", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsCheck(cons, scip->set, sol, checkintegrality, checklprows, printreason, result) );

   return SCIP_OKAY;
}

/** enforces single constraint for a given pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
SCIP_RETCODE SCIPenfopsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsAdded(cons));
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPenfopsCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnfops(cons, scip->set, solinfeasible, objinfeasible, result) );

   return SCIP_OKAY;
}

/** enforces single constraint for a given LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
SCIP_RETCODE SCIPenfolpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsAdded(cons));
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPenfolpCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnfolp(cons, scip->set, solinfeasible, result) );

   return SCIP_OKAY;
}

/** enforces single constraint for a given relaxation solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
SCIP_RETCODE SCIPenforelaxCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_SOL*             sol,                /**< solution to enforce */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsAdded(cons));
   assert(sol != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPenforelaxCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnforelax(cons, scip->set, sol, solinfeasible, result) );

   return SCIP_OKAY;
}

/** calls LP initialization method for single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
SCIP_RETCODE SCIPinitlpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to initialize */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility was detected while building the LP */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsAdded(cons));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinitlpCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsInitlp(cons, scip->set, infeasible) );

   return SCIP_OKAY;
}

/** calls separation method of single constraint for LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.
 */
SCIP_RETCODE SCIPsepalpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to separate */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsepalpCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSepalp(cons, scip->set, result) );

   return SCIP_OKAY;
}

/** calls separation method of single constraint for given primal solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.
 */
SCIP_RETCODE SCIPsepasolCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to separate */
   SCIP_SOL*             sol,                /**< primal solution that should be separated*/
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(sol != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsepasolCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSepasol(cons, scip->set, sol, result) );

   return SCIP_OKAY;
}

/** calls domain propagation method of single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.
 */
SCIP_RETCODE SCIPpropCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to propagate */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving loop */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPpropCons", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsProp(cons, scip->set, proptiming, result) );

   return SCIP_OKAY;
}

/** resolves propagation conflict of single constraint
 *
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
SCIP_RETCODE SCIPrespropCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to resolve conflict for */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information passed to the corresponding SCIPinferVarLbCons() or SCIPinferVarUbCons() call */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound which is sufficient to be explained */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsAdded(cons));
   assert(infervar != NULL);
   assert(bdchgidx != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPrespropCons", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsResprop(cons, scip->set, infervar, inferinfo, boundtype, bdchgidx, relaxedbd, result) );

   return SCIP_OKAY;
}

/** presolves of single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *
 *  @note This is an advanced method and should be used with caution.
 */
SCIP_RETCODE SCIPpresolCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to presolve */
   int                   nrounds,            /**< number of presolving rounds already done */
   SCIP_PRESOLTIMING     presoltiming,       /**< presolving timing(s) to be performed */
   int                   nnewfixedvars,      /**< number of variables fixed since the last call to the presolving method */
   int                   nnewaggrvars,       /**< number of variables aggregated since the last call to the presolving method */
   int                   nnewchgvartypes,    /**< number of variable type changes since the last call to the presolving method */
   int                   nnewchgbds,         /**< number of variable bounds tightened since the last call to the presolving method */
   int                   nnewholes,          /**< number of domain holes added since the last call to the presolving method */
   int                   nnewdelconss,       /**< number of deleted constraints since the last call to the presolving method */
   int                   nnewaddconss,       /**< number of added constraints since the last call to the presolving method */
   int                   nnewupgdconss,      /**< number of upgraded constraints since the last call to the presolving method */
   int                   nnewchgcoefs,       /**< number of changed coefficients since the last call to the presolving method */
   int                   nnewchgsides,       /**< number of changed left or right hand sides since the last call to the presolving method */
   int*                  nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to count total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to count total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to count total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgvartypes != NULL);
   assert(nchgbds != NULL);
   assert(naddholes != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPpresolCons", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsPresol(cons, scip->set, nrounds, presoltiming, nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
         nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars, nchgvartypes,
         nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides , result) );

   return SCIP_OKAY;
}

/** calls constraint activation notification method of single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *
 *  @note This is an advanced method and should be used with caution.  It may only be called for constraints that were not
 *      added to SCIP beforehand.
 */
SCIP_RETCODE SCIPactiveCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to notify */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsAdded(cons));
   assert(!SCIPconsIsDeleted(cons));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPactiveCons", FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsActive(cons, scip->set) );

   return SCIP_OKAY;
}

/** calls constraint deactivation notification method of single constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This is an advanced method and should be used with caution. It may only be called for constraints that were not
 *        added to SCIP beforehand.
 */
SCIP_RETCODE SCIPdeactiveCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to notify */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsAdded(cons));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPdeactiveCons", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsDeactive(cons, scip->set) );

   return SCIP_OKAY;
}

/** outputs constraint information to file stream via the message handler system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed.
 *  @note The file stream will not be flushed directly, this can be achieved by calling SCIPinfoMessage() printing a
 *        newline character.
 */
SCIP_RETCODE SCIPprintCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPconsPrint(cons, scip->set, scip->messagehdlr, file) );

   return SCIP_OKAY;
}

/** method to collect the variables of a constraint
 *
 *  If the number of variables is greater than the available slots in the variable array, nothing happens except that
 *  the success point is set to FALSE. With the method SCIPgetConsNVars() it is possible to get the number of variables
 *  a constraint has in its scope.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note The success pointer indicates if all variables were copied into the vars arrray.
 *
 *  @note It might be that a constraint handler does not support this functionality, in that case the success pointer is
 *        set to FALSE.
 */
SCIP_RETCODE SCIPgetConsVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the variables are wanted */
   SCIP_VAR**            vars,               /**< array to store the involved variable of the constraint */
   int                   varssize,           /**< available slots in vars array which is needed to check if the array is large enough */
   SCIP_Bool*            success             /**< pointer to store whether the variables are successfully copied */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetConsVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert(scip != NULL);
   assert(cons != NULL);
   assert(vars != NULL);
   assert(success != NULL);

   SCIP_CALL( SCIPconsGetVars(cons, scip->set, vars, varssize, success) );

   return SCIP_OKAY;
}

/** method to collect the number of variables of a constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note The success pointer indicates if the contraint handler was able to return the number of variables
 *
 *  @note It might be that a constraint handler does not support this functionality, in that case the success pointer is
 *        set to FALSE
 */
SCIP_RETCODE SCIPgetConsNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the number of variables is wanted */
   int*                  nvars,              /**< pointer to store the number of variables */
   SCIP_Bool*            success             /**< pointer to store whether the constraint successfully returned the number of variables */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetConsNVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nvars != NULL);
   assert(success != NULL);

   SCIP_CALL( SCIPconsGetNVars(cons, scip->set, nvars, success) );

   return SCIP_OKAY;
}
