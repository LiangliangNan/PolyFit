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
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objconshdlr.cpp
 * @brief  C++ wrapper for constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objconshdlr.h"




/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   scip::ObjConshdlr*    objconshdlr;        /**< constraint handler object */
   SCIP_Bool             deleteobject;       /**< should the constraint handler object be deleted when conshdlr is freed? */
};




/*
 * Callback methods of constraint handler
 */

extern "C"
{

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);
   assert(conshdlrdata->objconshdlr->scip_ != scip);

   if( conshdlrdata->objconshdlr->iscloneable() )
   {
      scip::ObjConshdlr* newobjconshdlr;
      newobjconshdlr = dynamic_cast<scip::ObjConshdlr*> (conshdlrdata->objconshdlr->clone(scip, valid));

      /* call include method of constraint handler object */
      SCIP_CALL( SCIPincludeObjConshdlr(scip, newobjconshdlr, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);
   assert(conshdlrdata->objconshdlr->scip_ == scip);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_free(scip, conshdlr) );

   /* free conshdlr object */
   if( conshdlrdata->deleteobject )
      delete conshdlrdata->objconshdlr;

   /* free conshdlr data */
   delete conshdlrdata;
   SCIPconshdlrSetData(conshdlr, NULL); /*lint !e64*/

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);
   assert(conshdlrdata->objconshdlr->scip_ == scip);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_init(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_exit(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_initpre(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_exitpre(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_initsol(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_exitsol(scip, conshdlr, conss, nconss, restart) );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_delete(scip, conshdlr, cons, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_trans(scip, conshdlr, sourcecons, targetcons) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_initlp(scip, conshdlr, conss, nconss, infeasible) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_sepalp(scip, conshdlr, conss, nconss, nusefulconss, result) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_sepasol(scip, conshdlr, conss, nconss, nusefulconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_enfolp(scip, conshdlr, conss, nconss, nusefulconss, solinfeasible, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_enforelax(scip, sol, conshdlr, conss, nconss, nusefulconss, solinfeasible, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_enfops(scip, conshdlr, conss, nconss, nusefulconss,
         solinfeasible, objinfeasible, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for primal solutions */
static
SCIP_DECL_CONSCHECK(consCheckObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_check(scip, conshdlr, conss, nconss, sol,
         checkintegrality, checklprows, printreason, completely, result) );

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_prop(scip, conshdlr, conss, nconss, nusefulconss, nmarkedconss, proptiming, result) );

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_presol(scip, conshdlr, conss, nconss, nrounds, presoltiming,
         nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
         nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
         nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
         ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_resprop(scip, conshdlr, cons, infervar, inferinfo, boundtype, bdchgidx,
         relaxedbd, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_lock(scip, conshdlr, cons, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_active(scip, conshdlr, cons) );

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_deactive(scip, conshdlr, cons) );

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_enable(scip, conshdlr, cons) );

   return SCIP_OKAY;
}


/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_disable(scip, conshdlr, cons) );

   return SCIP_OKAY;
}

/** variable deletion method of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelVarsObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_delvars(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_print(scip, conshdlr, cons, file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* sourceconshdlrdata;

   sourceconshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert(sourceconshdlrdata != NULL);
   assert(sourceconshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( sourceconshdlrdata->objconshdlr->scip_copy(scip, cons, name, sourcescip, sourceconshdlr, sourcecons, varmap, consmap,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_parse(scip, conshdlr, cons, name, str,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, success) );

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_getvars(scip, conshdlr, cons, vars, varssize, success) );

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_getnvars(scip, conshdlr, cons, nvars, success) );

   return SCIP_OKAY;
}

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsObj)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   SCIP_CALL( conshdlrdata->objconshdlr->scip_getdivebdchgs(scip, conshdlr, diveset, sol, success, infeasible) );

   return SCIP_OKAY;
}
}


/*
 * constraint handler specific interface methods
 */

/** creates the constraint handler for the given constraint handler object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjConshdlr*    objconshdlr,        /**< constraint handler object */
   SCIP_Bool             deleteobject        /**< should the constraint handler object be deleted when conshdlr is freed? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(objconshdlr != NULL);
   assert(objconshdlr->scip_ == scip);

   /* create obj constraint handler data */
   conshdlrdata = new SCIP_CONSHDLRDATA;
   conshdlrdata->objconshdlr = objconshdlr;
   conshdlrdata->deleteobject = deleteobject;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, objconshdlr->scip_name_, objconshdlr->scip_desc_,
         objconshdlr->scip_sepapriority_, objconshdlr->scip_enfopriority_, objconshdlr->scip_checkpriority_,
         objconshdlr->scip_sepafreq_, objconshdlr->scip_propfreq_, objconshdlr->scip_eagerfreq_,
         objconshdlr->scip_maxprerounds_,
         objconshdlr->scip_delaysepa_, objconshdlr->scip_delayprop_,
         objconshdlr->scip_needscons_, objconshdlr->scip_proptiming_, objconshdlr->scip_presoltiming_,
         conshdlrCopyObj,
         consFreeObj, consInitObj, consExitObj,
         consInitpreObj, consExitpreObj, consInitsolObj, consExitsolObj,
         consDeleteObj, consTransObj, consInitlpObj,
         consSepalpObj, consSepasolObj, consEnfolpObj, consEnforelaxObj, consEnfopsObj, consCheckObj,
         consPropObj, consPresolObj, consRespropObj, consLockObj,
         consActiveObj, consDeactiveObj,
         consEnableObj, consDisableObj, consDelVarsObj,
         consPrintObj, consCopyObj, consParseObj,
         consGetVarsObj, consGetNVarsObj, consGetDiveBdChgsObj, conshdlrdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the conshdlr object of the given name, or 0 if not existing */
scip::ObjConshdlr* SCIPfindObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, name);
   if( conshdlr == NULL )
      return 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->objconshdlr;
}

/** returns the conshdlr object for the given constraint handler */
scip::ObjConshdlr* SCIPgetObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->objconshdlr;
}
