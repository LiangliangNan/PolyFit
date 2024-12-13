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

/**@file   scip_prop.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for propagator plugins
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
#include "scip/prop.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_prop.h"
#include "scip/scip_prop.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates a propagator and includes it in SCIP.
 *
 *  @note method has all propagator callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludePropBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagator should be executed */
   int                   presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the propagator's presolving method */
   SCIP_DECL_PROPCOPY    ((*propcopy)),      /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   SCIP_DECL_PROPINITPRE ((*propinitpre)),   /**< presolving initialization method of propagator */
   SCIP_DECL_PROPEXITPRE ((*propexitpre)),   /**< presolving deinitialization method of propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_PROP* prop;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeProp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether propagator is already present */
   if( SCIPfindProp(scip, name) != NULL )
   {
      SCIPerrorMessage("propagator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpropCreate(&prop, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, freq, delay, timingmask, presolpriority, presolmaxrounds, presoltiming,
         propcopy, propfree, propinit, propexit, propinitpre, propexitpre, propinitsol, propexitsol,
         proppresol, propexec, propresprop, propdata) );
   SCIP_CALL( SCIPsetIncludeProp(scip->set, prop) );

   return SCIP_OKAY;
}

/** creates a propagator and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetPropInit(), SCIPsetPropExit(),
 *  SCIPsetPropCopy(), SCIPsetPropFree(), SCIPsetPropInitsol(), SCIPsetPropExitsol(),
 *  SCIPsetPropInitpre(), SCIPsetPropExitpre(), SCIPsetPropPresol(), and SCIPsetPropResprop().
 *
*  @note if you want to set all callbacks with a single method call, consider using SCIPincludeProp() instead
 */
SCIP_RETCODE SCIPincludePropBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP**           propptr,            /**< reference to a propagator pointer, or NULL */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagators should be executed */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_PROP* prop;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludePropBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether propagator is already present */
   if( SCIPfindProp(scip, name) != NULL )
   {
      SCIPerrorMessage("propagator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpropCreate(&prop, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, freq, delay, timingmask, 0, -1, SCIP_PRESOLTIMING_ALWAYS,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, propexec, NULL, propdata) );
   SCIP_CALL( SCIPsetIncludeProp(scip->set, prop) );

   if( propptr != NULL )
      *propptr = prop;

   return SCIP_OKAY;
}

/** sets copy method of propagator */
SCIP_RETCODE SCIPsetPropCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPCOPY    ((*propcopy))       /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetCopy(prop, propcopy);

   return SCIP_OKAY;
}

/** sets destructor method of propagator */
SCIP_RETCODE SCIPsetPropFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPFREE    ((*propfree))       /**< destructor of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetFree(prop, propfree);

   return SCIP_OKAY;
}

/** sets initialization method of propagator */
SCIP_RETCODE SCIPsetPropInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINIT    ((*propinit))       /**< initialize propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetInit(prop, propinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of propagator */
SCIP_RETCODE SCIPsetPropExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXIT    ((*propexit))       /**< deinitialize propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetExit(prop, propexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of propagator */
SCIP_RETCODE SCIPsetPropInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITSOL((*propinitsol))     /**< solving process initialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetInitsol(prop, propinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of propagator */
SCIP_RETCODE SCIPsetPropExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol))    /**< solving process deinitialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetExitsol(prop, propexitsol);

   return SCIP_OKAY;
}

/** sets preprocessing initialization method of propagator */
SCIP_RETCODE SCIPsetPropInitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITPRE((*propinitpre))     /**< preprocessing initialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropInitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetInitpre(prop, propinitpre);

   return SCIP_OKAY;
}

/** sets preprocessing deinitialization method of propagator */
SCIP_RETCODE SCIPsetPropExitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITPRE((*propexitpre))     /**< preprocessing deinitialization method of propagator */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropExitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetExitpre(prop, propexitpre);

   return SCIP_OKAY;
}

/** sets presolving method of propagator */
SCIP_RETCODE SCIPsetPropPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPPRESOL((*proppresol)),      /**< presolving method of propagator */
   int                   presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the propagator's presolving method */
   )
{
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropPresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);
   SCIP_CALL( SCIPpropSetPresol(prop, proppresol, presolpriority, presolmaxrounds, presoltiming) );

   name = SCIPpropGetName(prop);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, presolmaxrounds) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/presolpriority", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, presolpriority) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/presoltiming", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, (int) presoltiming) );

   return SCIP_OKAY;
}

/** sets propagation conflict resolving callback of propagator */
SCIP_RETCODE SCIPsetPropResprop(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop))    /**< propagation conflict resolving callback */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetPropResprop", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(prop != NULL);

   SCIPpropSetResprop(prop, propresprop);

   return SCIP_OKAY;
}


/** returns the propagator of the given name, or NULL if not existing */
SCIP_PROP* SCIPfindProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of propagator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindProp(scip->set, name);
}

/** returns the array of currently available propagators */
SCIP_PROP** SCIPgetProps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortProps(scip->set);

   return scip->set->props;
}

/** returns the number of currently available propagators */
int SCIPgetNProps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nprops;
}

/** sets the priority of a propagator */
SCIP_RETCODE SCIPsetPropPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   priority            /**< new priority of the propagator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPpropSetPriority(prop, scip->set, priority);

   return SCIP_OKAY;
}

/** sets the presolving priority of a propagator */
SCIP_RETCODE SCIPsetPropPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   presolpriority      /**< new presol priority of the propagator */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPpropSetPresolPriority(prop, scip->set, presolpriority);

   return SCIP_OKAY;
}
