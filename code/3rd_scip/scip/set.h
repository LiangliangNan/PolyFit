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

/**@file   set.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for global SCIP settings
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SET_H__
#define __SCIP_SET_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_bandit.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "scip/type_paramset.h"
#include "scip/type_event.h"
#include "scip/type_scip.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_cons.h"
#include "scip/type_disp.h"
#include "scip/type_heur.h"
#include "scip/type_compr.h"
#include "scip/type_nodesel.h"
#include "scip/type_presol.h"
#include "scip/type_pricer.h"
#include "scip/type_reader.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_table.h"
#include "scip/type_prop.h"
#include "scip/type_benders.h"
#include "scip/struct_set.h"


#ifdef NDEBUG
#include "scip/pub_misc.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** copies plugins from sourcescip to targetscip; in case that a constraint handler which does not need constraints
 *  cannot be copied, valid will return FALSE. All plugins can declare that, if their copy process failed, the 
 *  copied SCIP instance might not represent the same problem semantics as the original. 
 *  Note that in this case dual reductions might be invalid. */
SCIP_RETCODE SCIPsetCopyPlugins(
   SCIP_SET*             sourceset,          /**< source SCIP_SET data structure */
   SCIP_SET*             targetset,          /**< target SCIP_SET data structure */
   SCIP_Bool             copyreaders,        /**< should the file readers be copied */
   SCIP_Bool             copypricers,        /**< should the variable pricers be copied */
   SCIP_Bool             copyconshdlrs,      /**< should the constraint handlers be copied */
   SCIP_Bool             copyconflicthdlrs,  /**< should the conflict handlers be copied */
   SCIP_Bool             copypresolvers,     /**< should the presolvers be copied */
   SCIP_Bool             copyrelaxators,     /**< should the relaxators be copied */
   SCIP_Bool             copyseparators,     /**< should the separators be copied */
   SCIP_Bool             copycutselectors,   /**< should the cut selectors be copied */
   SCIP_Bool             copypropagators,    /**< should the propagators be copied */
   SCIP_Bool             copyheuristics,     /**< should the heuristics be copied */
   SCIP_Bool             copyeventhdlrs,     /**< should the event handlers be copied */
   SCIP_Bool             copynodeselectors,  /**< should the node selectors be copied */
   SCIP_Bool             copybranchrules,    /**< should the branchrules be copied */
   SCIP_Bool             copydisplays,       /**< should the display columns be copied */
   SCIP_Bool             copydialogs,        /**< should the dialogs be copied */
   SCIP_Bool             copytables,         /**< should the statistics tables be copied */
   SCIP_Bool             copyexprhdlrs,      /**< should the expression handlers be copied */
   SCIP_Bool             copynlpis,          /**< should the NLP interfaces be copied */
   SCIP_Bool*            allvalid            /**< pointer to store whether all plugins  were validly copied */
   );

/** copies parameters from sourcescip to targetscip */
SCIP_RETCODE SCIPsetCopyParams(
   SCIP_SET*             sourceset,          /**< source SCIP_SET data structure */
   SCIP_SET*             targetset,          /**< target SCIP_SET data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler of target SCIP */
   );

/** creates global SCIP settings */
SCIP_RETCODE SCIPsetCreate(
   SCIP_SET**            set,                /**< pointer to SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees global SCIP settings */
SCIP_RETCODE SCIPsetFree(
   SCIP_SET**            set,                /**< pointer to SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns current stage of SCIP */
SCIP_STAGE SCIPsetGetStage(
   SCIP_SET*             set                 /**< pointer to SCIP settings */
   );

/** creates a SCIP_Bool parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a SCIP_Longint parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a SCIP_Real parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** gets the fixing status value of an existing parameter */
SCIP_Bool SCIPsetIsParamFixed(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of the parameter */
   );

/** returns the pointer to the SCIP parameter with the given name */
SCIP_PARAM* SCIPsetGetParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of the parameter */
   );

/** gets the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetGetBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Int parameter */
SCIP_RETCODE SCIPsetGetIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   int*                  value               /**< pointer to store the parameter */
   );

/** gets the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPsetGetLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint*         value               /**< pointer to store the parameter */
   );

/** gets the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPsetGetRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Real*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Char parameter */
SCIP_RETCODE SCIPsetGetCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   char*                 value               /**< pointer to store the parameter */
   );

/** gets the value of an existing String parameter */
SCIP_RETCODE SCIPsetGetStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   char**                value               /**< pointer to store the parameter */
   );

/** changes the fixing status of an existing parameter */
SCIP_RETCODE SCIPsetChgParamFixed(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             fixed               /**< new fixing status of the parameter */
   );

/** changes the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetChgBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetSetBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   );

/** changes the default value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetSetDefaultBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             defaultvalue        /**< new default value of the parameter */
   );

/** changes the value of an existing Int parameter */
SCIP_RETCODE SCIPsetChgIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   int                   value               /**< new value of the parameter */
   );

/** changes the value of an existing Int parameter */
SCIP_RETCODE SCIPsetSetIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   int                   value               /**< new value of the parameter */
   );

/** changes the default value of an existing Int parameter */
SCIP_RETCODE SCIPsetSetDefaultIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   int                   defaultvalue        /**< new default value of the parameter */
   );

/** changes the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPsetChgLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPsetSetLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPsetChgRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPsetSetRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   SCIP_Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing Char parameter */
SCIP_RETCODE SCIPsetChgCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   char                  value               /**< new value of the parameter */
   );

/** changes the value of an existing Char parameter */
SCIP_RETCODE SCIPsetSetCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   char                  value               /**< new value of the parameter */
   );

/** changes the value of an existing String parameter */
SCIP_RETCODE SCIPsetChgStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   const char*           value               /**< new value of the parameter */
   );

/** changes the value of an existing String parameter */
SCIP_RETCODE SCIPsetSetStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   const char*           value               /**< new value of the parameter */
   );

/** reads parameters from a file */
SCIP_RETCODE SCIPsetReadParams(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< file name */
   );

/** writes all parameters in the parameter set to a file */
SCIP_RETCODE SCIPsetWriteParams(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   );

/** resets a single parameters to its default value */
SCIP_RETCODE SCIPsetResetParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name                /**< name of the parameter */
   );

/** resets all parameters to their default values */
SCIP_RETCODE SCIPsetResetParams(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** sets parameters to
 *
 *  - \ref SCIP_PARAMEMPHASIS_DEFAULT to use default values (see also SCIPsetResetParams())
 *  - \ref SCIP_PARAMEMPHASIS_COUNTER to get feasible and "fast" counting process
 *  - \ref SCIP_PARAMEMPHASIS_CPSOLVER to get CP like search (e.g. no LP relaxation)
 *  - \ref SCIP_PARAMEMPHASIS_EASYCIP to solve easy problems fast
 *  - \ref SCIP_PARAMEMPHASIS_FEASIBILITY to detect feasibility fast
 *  - \ref SCIP_PARAMEMPHASIS_HARDLP to be capable to handle hard LPs
 *  - \ref SCIP_PARAMEMPHASIS_OPTIMALITY to prove optimality fast
 *  - \ref SCIP_PARAMEMPHASIS_PHASEFEAS to find feasible solutions during a 3 phase solution process
 *  - \ref SCIP_PARAMEMPHASIS_PHASEIMPROVE to find improved solutions during a 3 phase solution process
 *  - \ref SCIP_PARAMEMPHASIS_PHASEPROOF to proof optimality during a 3 phase solution process
 */
SCIP_RETCODE SCIPsetSetEmphasis(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMEMPHASIS    paramemphasis,      /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** set parameters for reoptimization */
SCIP_RETCODE SCIPsetSetReoptimizationParams(
   SCIP_SET*             set,                /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** enable or disable all plugin timers depending on the value of the flag \p enabled */
void SCIPsetEnableOrDisablePluginClocks(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_Bool             enabled             /**< should plugin clocks be enabled? */
   );

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 */
SCIP_RETCODE SCIPsetSetSubscipsOff(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets heuristic parameters values to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for heuristic is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristic are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 */
SCIP_RETCODE SCIPsetSetHeuristics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets presolving parameters to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving is more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 */
SCIP_RETCODE SCIPsetSetPresolving(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets separating parameters to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the separating is done more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 */
SCIP_RETCODE SCIPsetSetSeparating(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** returns the array of all available SCIP parameters */
SCIP_PARAM** SCIPsetGetParams(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns the total number of all available SCIP parameters */
int SCIPsetGetNParams(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts file reader in file reader list */
SCIP_RETCODE SCIPsetIncludeReader(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_READER*          reader              /**< file reader */
   );

/** returns the file reader of the given name, or NULL if not existing */
SCIP_READER* SCIPsetFindReader(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of file reader */
   );

/** inserts variable pricer in variable pricer list */
SCIP_RETCODE SCIPsetIncludePricer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** returns the variable pricer of the given name, or NULL if not existing */
SCIP_PRICER* SCIPsetFindPricer(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of variable pricer */
   );

/** sorts pricers by priorities */
void SCIPsetSortPricers(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts pricers by name */
void SCIPsetSortPricersName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts Benders' decomposition into the Benders' decomposition list */
SCIP_RETCODE SCIPsetIncludeBenders(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** returns the Benders' decomposition of the given name, or NULL if not existing */
SCIP_BENDERS* SCIPsetFindBenders(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of Benders' decomposition */
   );

/** sorts Benders' decomposition by priorities */
void SCIPsetSortBenders(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts Benders' decomposition by name */
void SCIPsetSortBendersName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts constraint handler in constraint handler list */
SCIP_RETCODE SCIPsetIncludeConshdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** reinserts a constraint handler with modified sepa priority into the sepa priority sorted array */
void SCIPsetReinsertConshdlrSepaPrio(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler to be reinserted */
   int                   oldpriority         /**< the old separation priority of constraint handler */
   );

/** returns the constraint handler of the given name, or NULL if not existing */
SCIP_CONSHDLR* SCIPsetFindConshdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of constraint handler */
   );

/** inserts conflict handler in conflict handler list */
SCIP_RETCODE SCIPsetIncludeConflicthdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** returns the conflict handler of the given name, or NULL if not existing */
SCIP_CONFLICTHDLR* SCIPsetFindConflicthdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of conflict handler */
   );

/** sorts conflict handlers by priorities */
void SCIPsetSortConflicthdlrs(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts conflict handlers by name */
void SCIPsetSortConflicthdlrsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts presolver in presolver list */
SCIP_RETCODE SCIPsetIncludePresol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** returns the presolver of the given name, or NULL if not existing */
SCIP_PRESOL* SCIPsetFindPresol(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of presolver */
   );

/** sorts presolvers by priorities */
void SCIPsetSortPresols(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts presolvers by name */
void SCIPsetSortPresolsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts relaxator in relaxator list */
SCIP_RETCODE SCIPsetIncludeRelax(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAX*           relax               /**< relaxator */
   );

/** returns the relaxator of the given name, or NULL if not existing */
SCIP_RELAX* SCIPsetFindRelax(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of relaxator */
   );

/** sorts relaxators by priorities */
void SCIPsetSortRelaxs(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts relaxators by name */
void SCIPsetSortRelaxsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts separator in separator list */
SCIP_RETCODE SCIPsetIncludeSepa(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SEPA*            sepa                /**< separator */
   );

/** returns the separator of the given name, or NULL if not existing */
SCIP_SEPA* SCIPsetFindSepa(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of separator */
   );

/** sorts separators by priorities */
void SCIPsetSortSepas(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts separators by name */
void SCIPsetSortSepasName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts cut selector in cut selector list */
SCIP_RETCODE SCIPsetIncludeCutsel(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** returns the cut selector of the given name, or NULL if not existing */
SCIP_CUTSEL* SCIPsetFindCutsel(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of separator */
   );

/** sorts cut selectors by priorities */
void SCIPsetSortCutsels(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts propagator in propagator list */
SCIP_RETCODE SCIPsetIncludeProp(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROP*            prop                /**< propagator */
   );

/** returns the propagator of the given name, or NULL if not existing */
SCIP_PROP* SCIPsetFindProp(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of propagator */
   );

/** sorts propagators by priorities */
void SCIPsetSortProps(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts propagators by priorities for presolving */
void SCIPsetSortPropsPresol(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts propagators w.r.t. names */
void SCIPsetSortPropsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts concurrent solver type into the concurrent solver type list */
SCIP_RETCODE SCIPsetIncludeConcsolverType(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   );

/** returns the concurrent solver type with the given name, or NULL if not existing */
SCIP_CONCSOLVERTYPE* SCIPsetFindConcsolverType(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of concurrent solver type */
   );

/** inserts concurrent solver into the concurrent solver list */
SCIP_RETCODE SCIPsetIncludeConcsolver(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** frees all concurrent solvers in the concurrent solver list */
SCIP_RETCODE SCIPsetFreeConcsolvers(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts primal heuristic in primal heuristic list */
SCIP_RETCODE SCIPsetIncludeHeur(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** returns the primal heuristic of the given name, or NULL if not existing */
SCIP_HEUR* SCIPsetFindHeur(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of primal heuristic */
   );

/** sorts heuristics by priorities */
void SCIPsetSortHeurs(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts heuristics by name */
void SCIPsetSortHeursName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts tree compression in tree compression list */
SCIP_RETCODE SCIPsetIncludeCompr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COMPR*           compr               /**< tree compression */
   );

/** returns the tree compression of the given name, or NULL if not existing */
SCIP_COMPR* SCIPsetFindCompr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of tree compression */
   );

/** sorts compressions by priorities */
void SCIPsetSortComprs(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts heuristics by names */
void SCIPsetSortComprsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts event handler in event handler list */
SCIP_RETCODE SCIPsetIncludeEventhdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** returns the event handler of the given name, or NULL if not existing */
SCIP_EVENTHDLR* SCIPsetFindEventhdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of event handler */
   );

/** inserts node selector in node selector list */
SCIP_RETCODE SCIPsetIncludeNodesel(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** returns the node selector of the given name, or NULL if not existing */
SCIP_NODESEL* SCIPsetFindNodesel(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of event handler */
   );

/** returns node selector with highest priority in the current mode */
SCIP_NODESEL* SCIPsetGetNodesel(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** inserts branching rule in branching rule list */
SCIP_RETCODE SCIPsetIncludeBranchrule(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** returns the branching rule of the given name, or NULL if not existing */
SCIP_BRANCHRULE* SCIPsetFindBranchrule(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of event handler */
   );

/** sorts branching rules by priorities */
void SCIPsetSortBranchrules(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sorts branching rules by name */
void SCIPsetSortBranchrulesName(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts display column in display column list */
SCIP_RETCODE SCIPsetIncludeDisp(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DISP*            disp                /**< display column */
   );

/** returns the display column of the given name, or NULL if not existing */
SCIP_DISP* SCIPsetFindDisp(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of display */
   );

/** inserts statistics table in statistics table list */
SCIP_RETCODE SCIPsetIncludeTable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TABLE*           table               /**< statistics table */
   );

/** returns the statistics table of the given name, or NULL if not existing */
SCIP_TABLE* SCIPsetFindTable(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of statistics table */
   );

/** inserts dialog in dialog list */
SCIP_RETCODE SCIPsetIncludeDialog(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** returns if the dialog already exists */
SCIP_Bool SCIPsetExistsDialog(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** inserts expression handler in expression handler list */
SCIP_RETCODE SCIPsetIncludeExprhdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns the expression handler of the given name, or NULL if not existing */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPfindExprhdlr() macro */
SCIP_EXPRHDLR* SCIPsetFindExprhdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of expression handler */
   );

/** sorts expression handlers by name */
void SCIPsetSortExprhdlrs(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** inserts NLPI in NLPI list */
SCIP_RETCODE SCIPsetIncludeNlpi(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi                /**< NLPI */
   );

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_NLPI* SCIPsetFindNlpi(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of NLPI */
   );

/** sorts NLPIs by priorities */
void SCIPsetSortNlpis(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** set priority of an NLPI */
void SCIPsetSetPriorityNlpi(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of NLPI */
   );

/** inserts information about an external code in external codes list */
SCIP_RETCODE SCIPsetIncludeExternalCode(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of external code */
   const char*           description         /**< description of external code, can be NULL */
   );

/** inserts bandit virtual function table into set */
SCIP_RETCODE SCIPsetIncludeBanditvtable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BANDITVTABLE*    banditvtable        /**< bandit algorithm virtual function table */
   );

/** returns the bandit virtual function table of the given name, or NULL if not existing */
SCIP_BANDITVTABLE* SCIPsetFindBanditvtable(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of bandit algorithm virtual function table */
   );

/** calls init methods of all plugins */
SCIP_RETCODE SCIPsetInitPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** calls exit methods of all plugins */
SCIP_RETCODE SCIPsetExitPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** calls initpre methods of all plugins */
SCIP_RETCODE SCIPsetInitprePlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** calls exitpre methods of all plugins */
SCIP_RETCODE SCIPsetExitprePlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** calls initsol methods of all plugins */
SCIP_RETCODE SCIPsetInitsolPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** calls exitsol methods of all plugins */
SCIP_RETCODE SCIPsetExitsolPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   );

/** calculate memory size for dynamically allocated arrays */
int SCIPsetCalcMemGrowSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** calculate memory size for tree array */
int SCIPsetCalcTreeGrowSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** calculate memory size for path array */
int SCIPsetCalcPathGrowSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** sets verbosity level for message output */
SCIP_RETCODE SCIPsetSetVerbLevel(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VERBLEVEL        verblevel           /**< verbosity level for message output */
   );

/** sets feasibility tolerance */
SCIP_RETCODE SCIPsetSetFeastol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< LP data, or NULL */
   SCIP_Real             feastol             /**< new feasibility tolerance */
   );

/** sets feasibility tolerance for reduced costs in LP solution */
SCIP_RETCODE SCIPsetSetDualfeastol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             dualfeastol         /**< new reduced costs feasibility tolerance */
   );

/** sets LP convergence tolerance used in barrier algorithm */
SCIP_RETCODE SCIPsetSetBarrierconvtol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   );

/** sets primal feasibility tolerance for relaxations (relaxfeastol)
 *
 * @note Set to SCIP_INVALID to apply relaxation-specific feasibility tolerance only.
 *
 * @return Previous value of relaxfeastol.
 */
SCIP_Real SCIPsetSetRelaxfeastol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             relaxfeastol        /**< new primal feasibility tolerance for relaxations, or SCIP_INVALID */
   );

/** marks that some limit parameter was changed */
void SCIPsetSetLimitChanged(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns the maximal number of variables priced into the LP per round */
int SCIPsetGetPriceMaxvars(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             root                /**< are we at the root node? */
   );

/** returns the maximal number of cuts separated per round */
int SCIPsetGetSepaMaxcuts(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             root                /**< are we at the root node? */
   );

/** returns the maximal ratio between coefficients to ensure in rowprep cleanup */
SCIP_Real SCIPsetGetSepaMaxCoefRatioRowprep(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns user defined objective value (in original space) for reference purposes */
SCIP_Real SCIPsetGetReferencevalue(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns debug solution data */
SCIP_DEBUGSOLDATA* SCIPsetGetDebugSolData(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** Checks, if an iteratively updated value is reliable or should be recomputed from scratch.
 *  This is useful, if the value, e.g., the activity of a linear constraint or the pseudo objective value, gets a high
 *  absolute value during the optimization process which is later reduced significantly. In this case, the last digits
 *  were cancelled out when increasing the value and are random after decreasing it.
 *  We dot not consider the cancellations which can occur during increasing the absolute value because they just cannot
 *  be expressed using fixed precision floating point arithmetic, anymore.
 *  The idea to get more reliable values is to always store the last reliable value, where increasing the absolute of
 *  the value is viewed as preserving reliability. Then, after each update, the new absolute value can be compared
 *  against the last reliable one with this method, checking whether it was decreased by a factor of at least
 *  "lp/recompfac" and should be recomputed.
 */
SCIP_Bool SCIPsetIsUpdateUnreliable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newvalue,           /**< new value after update */
   SCIP_Real             oldvalue            /**< old value, i.e., last reliable value */
   );

/** modifies an initial seed value with the global shift of random seeds */
unsigned int SCIPsetInitializeRandomSeed(
   SCIP_SET*             set,                /**< global SCIP settings */
   unsigned int          initialseedvalue    /**< initial seed value to be modified */
   );

/** returns value treated as infinity */
SCIP_Real SCIPsetInfinity(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns the minimum value that is regarded as huge and should be handled separately (e.g., in activity
 *  computation)
 */
SCIP_Real SCIPsetGetHugeValue(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns value treated as zero */
SCIP_Real SCIPsetEpsilon(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns value treated as zero for sums of floating point values */
SCIP_Real SCIPsetSumepsilon(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns feasibility tolerance for constraints */
#ifdef __GNUC__
__attribute__ ((pure))
#endif
SCIP_Real SCIPsetFeastol(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns factor w.r.t. primal feasibility tolerance that determines default (and maximal) feasibility tolerance */
SCIP_Real SCIPsetLPFeastolFactor(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns feasibility tolerance for reduced costs */
#ifdef __GNUC__
__attribute__ ((pure))
#endif
SCIP_Real SCIPsetDualfeastol(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns convergence tolerance used in barrier algorithm */
SCIP_Real SCIPsetBarrierconvtol(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns minimal variable distance value to use for pseudo cost updates */
SCIP_Real SCIPsetPseudocosteps(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns minimal minimal objective distance value to use for pseudo cost updates */
SCIP_Real SCIPsetPseudocostdelta(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** return the delta to use for computing the cutoff bound for integral objectives */
SCIP_Real SCIPsetCutoffbounddelta(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** return the primal feasibility tolerance for relaxations */
SCIP_Real SCIPsetRelaxfeastol(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns minimal decrease factor that causes the recomputation of a value
 *  (e.g., pseudo objective) instead of an update */
SCIP_Real SCIPsetRecompfac(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** checks, if values are in range of epsilon */
SCIP_Bool SCIPsetIsEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) lower than val2 */
SCIP_Bool SCIPsetIsLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) greater than val2 */
SCIP_Bool SCIPsetIsLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) greater than val2 */
SCIP_Bool SCIPsetIsGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) lower than val2 */
SCIP_Bool SCIPsetIsGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is (positive) infinite */
SCIP_Bool SCIPsetIsInfinity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against infinity */
   );

/** checks, if value is huge and should be handled separately (e.g., in activity computation) */
SCIP_Bool SCIPsetIsHugeValue(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be checked whether it is huge */
   );

/** checks, if value is in range epsilon of 0.0 */
SCIP_Bool SCIPsetIsZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than epsilon */
SCIP_Bool SCIPsetIsPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -epsilon */
SCIP_Bool SCIPsetIsNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is integral within epsilon */
SCIP_Bool SCIPsetIsIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
SCIP_Bool SCIPsetIsScalingIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val,                /**< unscaled value to check for scaled integrality */
   SCIP_Real             scalar              /**< value to scale val with for checking for integrality */
   );

/** checks, if given fractional part is smaller than epsilon */
SCIP_Bool SCIPsetIsFracIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value + feasibility tolerance down to the next integer in epsilon tolerance */
SCIP_Real SCIPsetFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value - feasibility tolerance up to the next integer in epsilon tolerance */
SCIP_Real SCIPsetCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value to the nearest integer in epsilon tolerance */
SCIP_Real SCIPsetRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
SCIP_Real SCIPsetFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to return fractional part for */
   );

/** checks, if values are in range of sumepsilon */
SCIP_Bool SCIPsetIsSumEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPsetIsSumLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPsetIsSumLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPsetIsSumGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPsetIsSumGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range sumepsilon of 0.0 */
SCIP_Bool SCIPsetIsSumZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than sumepsilon */
SCIP_Bool SCIPsetIsSumPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -sumepsilon */
SCIP_Bool SCIPsetIsSumNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value + sumepsilon tolerance down to the next integer */
SCIP_Real SCIPsetSumFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value - sumepsilon tolerance up to the next integer */
SCIP_Real SCIPsetSumCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value to the nearest integer in sumepsilon tolerance */
SCIP_Real SCIPsetSumRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   );

/** returns fractional part of value, i.e. x - floor(x) in sumepsilon tolerance */
SCIP_Real SCIPsetSumFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if relative difference of values is in range of feastol */
SCIP_Bool SCIPsetIsFeasEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than feastol */
SCIP_Bool SCIPsetIsFeasLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than feastol */
SCIP_Bool SCIPsetIsFeasLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than feastol */
SCIP_Bool SCIPsetIsFeasGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
SCIP_Bool SCIPsetIsFeasGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range feasibility tolerance of 0.0 */
SCIP_Bool SCIPsetIsFeasZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than feasibility tolerance */
SCIP_Bool SCIPsetIsFeasPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -feasibility tolerance */
SCIP_Bool SCIPsetIsFeasNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is integral within the feasibility bounds */
SCIP_Bool SCIPsetIsFeasIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if given fractional part is smaller than feastol */
SCIP_Bool SCIPsetIsFeasFracIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value + feasibility tolerance down to the next integer */
SCIP_Real SCIPsetFeasFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value - feasibility tolerance up to the next integer */
SCIP_Real SCIPsetFeasCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value to the nearest integer in feasibility tolerance */
SCIP_Real SCIPsetFeasRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. x - floor(x) in feasibility tolerance */
SCIP_Real SCIPsetFeasFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to return fractional part for */
   );

/** checks, if relative difference of values is in range of dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range dual feasibility tolerance of 0.0 */
SCIP_Bool SCIPsetIsDualfeasZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is integral within the dual feasibility bounds */
SCIP_Bool SCIPsetIsDualfeasIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if given fractional part is smaller than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasFracIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value + dual feasibility tolerance down to the next integer */
SCIP_Real SCIPsetDualfeasFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value - dual feasibility tolerance up to the next integer */
SCIP_Real SCIPsetDualfeasCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** rounds value to the nearest integer in dual feasibility tolerance */
SCIP_Real SCIPsetDualfeasRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. x - floor(x) in dual feasibility tolerance */
SCIP_Real SCIPsetDualfeasFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to return fractional part for */
   );

/** checks, if the given new lower bound is at least min(oldub - oldlb, |oldlb|) times the bound
 *  strengthening epsilon better than the old one or the change in the lower bound would fix the
 *  sign of the variable
 */
SCIP_Bool SCIPsetIsLbBetter(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   );

/** checks, if the given new upper bound is at least min(oldub - oldlb, |oldub|) times the bound
 *  strengthening epsilon better than the old one or the change in the upper bound would fix the
 *  sign of the variable
 */
SCIP_Bool SCIPsetIsUbBetter(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   );

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy */
SCIP_Bool SCIPsetIsEfficacious(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             root,               /**< should the root's minimal cut efficacy be used? */
   SCIP_Real             efficacy            /**< efficacy of the cut */
   );

/** checks, if relative difference of values is in range of epsilon */
SCIP_Bool SCIPsetIsRelEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than epsilon */
SCIP_Bool SCIPsetIsRelLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
SCIP_Bool SCIPsetIsRelLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than epsilon */
SCIP_Bool SCIPsetIsRelGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
SCIP_Bool SCIPsetIsRelGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of values is in range of sumepsilon */
SCIP_Bool SCIPsetIsSumRelEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
SCIP_Bool SCIPsetIsSumRelLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
SCIP_Bool SCIPsetIsSumRelLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
SCIP_Bool SCIPsetIsSumRelGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
SCIP_Bool SCIPsetIsSumRelGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** returns the flag indicating whether sub-SCIPs that could cause recursion have been deactivated */
SCIP_Bool SCIPsetGetSubscipsOff(
   SCIP_SET*             set                 /**< global SCIP settings */
   );


#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsetInfinity(set)               ( (set)->num_infinity )
#define SCIPsetGetHugeValue(set)           ( (set)->num_hugeval )
#define SCIPsetEpsilon(set)                ( (set)->num_epsilon )
#define SCIPsetSumepsilon(set)             ( (set)->num_sumepsilon )
#define SCIPsetFeastol(set)                ( (set)->num_feastol )
#define SCIPsetLPFeastolFactor(set)        ( (set)->num_lpfeastolfactor )
#define SCIPsetDualfeastol(set)            ( (set)->num_dualfeastol )
#define SCIPsetBarrierconvtol(set)         ( (set)->num_barrierconvtol )
#define SCIPsetPseudocosteps(set)          ( (set)->num_pseudocosteps )
#define SCIPsetPseudocostdelta(set)        ( (set)->num_pseudocostdelta )
#define SCIPsetRelaxfeastol(set)           ( (set)->num_relaxfeastol )
#define SCIPsetCutoffbounddelta(set)       ( MIN(100.0 * SCIPsetFeastol(set), 0.0001) )
#define SCIPsetRecompfac(set)              ( (set)->num_recompfac )
#define SCIPsetIsEQ(set, val1, val2)       ( EPSEQ(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsLT(set, val1, val2)       ( EPSLT(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsLE(set, val1, val2)       ( EPSLE(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsGT(set, val1, val2)       ( EPSGT(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsGE(set, val1, val2)       ( EPSGE(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsInfinity(set, val)        ( (val) >= (set)->num_infinity )
#define SCIPsetIsHugeValue(set, val)       ( (val) >= (set)->num_hugeval )
#define SCIPsetIsZero(set, val)            ( EPSZ(val, (set)->num_epsilon) )
#define SCIPsetIsPositive(set, val)        ( EPSP(val, (set)->num_epsilon) )
#define SCIPsetIsNegative(set, val)        ( EPSN(val, (set)->num_epsilon) )
#define SCIPsetIsIntegral(set, val)        ( EPSISINT(val, (set)->num_epsilon) )
#define SCIPsetIsScalingIntegral(set, val, scalar)                      \
   ( EPSISINT((scalar)*(val), MAX(REALABS(scalar), 1.0)*(set)->num_epsilon) )
#define SCIPsetIsFracIntegral(set, val)    ( !EPSP(val, (set)->num_epsilon) )
#define SCIPsetFloor(set, val)             ( EPSFLOOR(val, (set)->num_epsilon) )
#define SCIPsetCeil(set, val)              ( EPSCEIL(val, (set)->num_epsilon) )
#define SCIPsetRound(set, val)             ( EPSROUND(val, (set)->num_epsilon) )
#define SCIPsetFrac(set, val)              ( EPSFRAC(val, (set)->num_epsilon) )

#define SCIPsetIsSumEQ(set, val1, val2)    ( EPSEQ(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumLT(set, val1, val2)    ( EPSLT(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumLE(set, val1, val2)    ( EPSLE(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumGT(set, val1, val2)    ( EPSGT(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumGE(set, val1, val2)    ( EPSGE(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumZero(set, val)         ( EPSZ(val, (set)->num_sumepsilon) )
#define SCIPsetIsSumPositive(set, val)     ( EPSP(val, (set)->num_sumepsilon) )
#define SCIPsetIsSumNegative(set, val)     ( EPSN(val, (set)->num_sumepsilon) )
#define SCIPsetSumFloor(set, val)          ( EPSFLOOR(val, (set)->num_sumepsilon) )
#define SCIPsetSumCeil(set, val)           ( EPSCEIL(val, (set)->num_sumepsilon) )
#define SCIPsetSumRound(set, val)          ( EPSROUND(val, (set)->num_sumepsilon) )
#define SCIPsetSumFrac(set, val)           ( EPSFRAC(val, (set)->num_sumepsilon) )

#define SCIPsetIsFeasEQ(set, val1, val2)   ( EPSZ(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasLT(set, val1, val2)   ( EPSN(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasLE(set, val1, val2)   ( !EPSP(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasGT(set, val1, val2)   ( EPSP(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasGE(set, val1, val2)   ( !EPSN(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasZero(set, val)        ( EPSZ(val, (set)->num_feastol) )
#define SCIPsetIsFeasPositive(set, val)    ( EPSP(val, (set)->num_feastol) )
#define SCIPsetIsFeasNegative(set, val)    ( EPSN(val, (set)->num_feastol) )
#define SCIPsetIsFeasIntegral(set, val)    ( EPSISINT(val, (set)->num_feastol) )
#define SCIPsetIsFeasFracIntegral(set, val) ( !EPSP(val, (set)->num_feastol) )
#define SCIPsetFeasFloor(set, val)         ( EPSFLOOR(val, (set)->num_feastol) )
#define SCIPsetFeasCeil(set, val)          ( EPSCEIL(val, (set)->num_feastol) )
#define SCIPsetFeasRound(set, val)         ( EPSROUND(val, (set)->num_feastol) )
#define SCIPsetFeasFrac(set, val)          ( EPSFRAC(val, (set)->num_feastol) )

#define SCIPsetIsDualfeasEQ(set, val1, val2)   ( EPSZ(SCIPrelDiff(val1, val2), (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasLT(set, val1, val2)   ( EPSN(SCIPrelDiff(val1, val2), (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasLE(set, val1, val2)   ( !EPSP(SCIPrelDiff(val1, val2), (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasGT(set, val1, val2)   ( EPSP(SCIPrelDiff(val1, val2), (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasGE(set, val1, val2)   ( !EPSN(SCIPrelDiff(val1, val2), (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasZero(set, val)        ( EPSZ(val, (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasPositive(set, val)    ( EPSP(val, (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasNegative(set, val)    ( EPSN(val, (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasIntegral(set, val)    ( EPSISINT(val, (set)->num_dualfeastol) )
#define SCIPsetIsDualfeasFracIntegral(set, val) ( !EPSP(val, (set)->num_dualfeastol) )
#define SCIPsetDualfeasFloor(set, val)         ( EPSFLOOR(val, (set)->num_dualfeastol) )
#define SCIPsetDualfeasCeil(set, val)          ( EPSCEIL(val, (set)->num_dualfeastol) )
#define SCIPsetDualfeasRound(set, val)         ( EPSROUND(val, (set)->num_dualfeastol) )
#define SCIPsetDualfeasFrac(set, val)          ( EPSFRAC(val, (set)->num_dualfeastol) )

#define SCIPsetIsLbBetter(set, newlb, oldlb, oldub) ( ((oldlb) < 0.0 && (newlb) >= 0.0) || EPSGT(newlb, oldlb, \
         set->num_boundstreps * MAX(MIN((oldub) - (oldlb), REALABS(oldlb)), 1e-3)) )
#define SCIPsetIsUbBetter(set, newub, oldlb, oldub) ( ((oldub) > 0.0 && (newub) <= 0.0) || EPSLT(newub, oldub, \
         set->num_boundstreps * MAX(MIN((oldub) - (oldlb), REALABS(oldub)), 1e-3)) )
#define SCIPsetIsEfficacious(set, root, efficacy) \
   ( root ? EPSP(efficacy, (set)->sepa_minefficacyroot) : EPSP(efficacy, (set)->sepa_minefficacy) )

#define SCIPsetIsRelEQ(set, val1, val2)    ( EPSZ(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelLT(set, val1, val2)    ( EPSN(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelLE(set, val1, val2)    ( !EPSP(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelGT(set, val1, val2)    ( EPSP(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelGE(set, val1, val2)    ( !EPSN(SCIPrelDiff(val1, val2), (set)->num_epsilon) )

#define SCIPsetIsSumRelEQ(set, val1, val2) ( EPSZ(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelLT(set, val1, val2) ( EPSN(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelLE(set, val1, val2) ( !EPSP(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelGT(set, val1, val2) ( EPSP(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelGE(set, val1, val2) ( !EPSN(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsUpdateUnreliable(set, newvalue, oldvalue) \
   ( (ABS(oldvalue) / MAX(ABS(newvalue), set->num_epsilon)) >= set->num_recompfac )
#define SCIPsetInitializeRandomSeed(set, val) ( (val + (set)->random_randomseedshift) )

#define SCIPsetGetSubscipsOff(set)         ( (set)->subscipsoff )

#endif

#define SCIPsetAllocBuffer(set,ptr)             ( (BMSallocBufferMemory((set)->buffer, (ptr)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetAllocBufferSize(set,ptr,size)    ( (BMSallocBufferMemorySize((set)->buffer, (ptr), (size)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetAllocBufferArray(set,ptr,num)    ( (BMSallocBufferMemoryArray((set)->buffer, (ptr), (num)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetDuplicateBufferSize(set,ptr,source,size) ( (BMSduplicateBufferMemory((set)->buffer, (ptr), (source), (size)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetDuplicateBufferArray(set,ptr,source,num) ( (BMSduplicateBufferMemoryArray((set)->buffer, (ptr), (source), (num)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetReallocBufferSize(set,ptr,size)  ( (BMSreallocBufferMemorySize((set)->buffer, (ptr), (size)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetReallocBufferArray(set,ptr,num)  ( (BMSreallocBufferMemoryArray((set)->buffer, (ptr), (num)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetFreeBuffer(set,ptr)              BMSfreeBufferMemory((set)->buffer, (ptr))
#define SCIPsetFreeBufferSize(set,ptr)          BMSfreeBufferMemorySize((set)->buffer, (ptr))
#define SCIPsetFreeBufferArray(set,ptr)         BMSfreeBufferMemoryArray((set)->buffer, (ptr))

#define SCIPsetAllocCleanBuffer(set,ptr)             ( (BMSallocBufferMemory((set)->cleanbuffer, (ptr)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetAllocCleanBufferSize(set,ptr,size)    ( (BMSallocBufferMemorySize((set)->cleanbuffer, (ptr), (size)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetAllocCleanBufferArray(set,ptr,num)    ( (BMSallocBufferMemoryArray((set)->cleanbuffer, (ptr), (num)) == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPsetFreeCleanBuffer(set,ptr)              BMSfreeBufferMemory((set)->cleanbuffer, (ptr))
#define SCIPsetFreeCleanBufferSize(set,ptr)          BMSfreeBufferMemorySize((set)->cleanbuffer, (ptr))
#define SCIPsetFreeCleanBufferArray(set,ptr)         BMSfreeBufferMemoryArray((set)->cleanbuffer, (ptr))

/* if we have a C99 compiler */
#ifdef SCIP_HAVE_VARIADIC_MACROS

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPsetDebugMsg(set, ...)       SCIPsetPrintDebugMessage(set, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPsetDebugMsgPrint(set, ...)  SCIPsetDebugMessagePrint(set, __VA_ARGS__)
#else
#define SCIPsetDebugMsg(set, ...)       while ( FALSE ) SCIPsetPrintDebugMessage(set, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPsetDebugMsgPrint(set, ...)  while ( FALSE ) SCIPsetDebugMessagePrint(set, __VA_ARGS__)
#endif

#else
/* if we do not have a C99 compiler, use a workaround that prints a message, but not the file and linenumber */

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPsetDebugMsg                 printf("debug: "), SCIPsetDebugMessagePrint
#define SCIPsetDebugMsgPrint            printf("debug: "), SCIPsetDebugMessagePrint
#else
#define SCIPsetDebugMsg                 while ( FALSE ) SCIPsetDebugMsgPrint
#define SCIPsetDebugMsgPrint            while ( FALSE ) SCIPsetDebugMessagePrint
#endif

#endif


/** prints a debug message */
#ifdef __GNUC__
__attribute__((format(printf, 4, 5)))
#endif
SCIP_EXPORT
void SCIPsetPrintDebugMessage(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a debug message without precode */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT
void SCIPsetDebugMessagePrint(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );


#ifdef __cplusplus
}
#endif

#endif
