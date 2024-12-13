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

/**@file   dialog_default.h
 * @ingroup DIALOGS
 * @brief  default user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DIALOG_DEFAULT_H__
#define __SCIP_DIALOG_DEFAULT_H__

#include "scip/def.h"
#include "scip/type_dialog.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup DIALOGS
 *
 * @{
 */

/** standard menu dialog execution method, that displays it's help screen if the remaining command line is empty */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenu);

/** standard menu dialog execution method, that doesn't display it's help screen */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenuLazy);

/** dialog execution method for the change add constraint */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeAddCons);

/** dialog execution method for the change bounds command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeBounds);

/** dialog execution method for the freetransproblem command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeFreetransproblem);

/** dialog execution method for the changing the objective sense */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeObjSense);

/** dialog execution method for the checksol command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChecksol);

/** dialog execution method for the cliquegraph command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCliquegraph);

/** dialog execution method for the display benders command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayBenders);

/** dialog execution method for the display branching command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayBranching);

/** dialog execution method for the display compression command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayCompression);

/** dialog execution method for the display conflict command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConflict);

/** dialog execution method for the display conshdlrs command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConshdlrs);

/** dialog execution method for the display displaycols command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayDisplaycols);

/** dialog execution method for the display exprhdlrs command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayExprhdlrs);

/** dialog execution method for the display cutselectors command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayCutselectors);

/** dialog execution method for the display heuristics command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayHeuristics);

/** dialog execution method for the display memory command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayMemory);

/** dialog execution method for the display nodeselectors command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNodeselectors);

/** dialog execution method for the display nlpi command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNlpi);

/** dialog execution method for the display parameters command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayParameters);

/** dialog execution method for the display presolvers command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPresolvers);

/** dialog execution method for the display pricer command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPricers);

/** dialog execution method for the display problem command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayProblem);

/** dialog execution method for the display propagators command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPropagators);

/** dialog execution method for the display readers command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayReaders);

/** dialog execution method for the display relaxators command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayRelaxators);

/** dialog execution method for the display separators command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySeparators);

/** dialog execution method for the display solution command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySolution);

/** dialog execution method for the display finitesolution command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayFiniteSolution);

/** dialog execution method for the display dual solution command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayDualSolution);

/** dialog execution method for the display of solutions in the pool command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySolutionPool);

/** dialog execution method for the display subproblem command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySubproblem);

/** dialog execution method for the display subsolution command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySubSolution);

/** dialog execution method for the display statistics command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayStatistics);

/** dialog execution method for the display reoptstatistics command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayReoptStatistics);

/** dialog execution method for the display transproblem command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTransproblem);

/** dialog execution method for the display value command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayValue);

/** dialog execution method for the display varbranchstatistics command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayVarbranchstatistics);

/** dialog execution method for the display LP solution quality command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayLPSolutionQuality);

/** dialog execution method for the display transsolution command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTranssolution);

/** dialog execution method for the help command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecHelp);

/** dialog execution method for the free command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecFree);

/** dialog execution method for the newstart command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecNewstart);

/** dialog execution method for the transform command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecTransform);

/** dialog execution method for the optimize command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecOptimize);

/** dialog execution method for the parallelopt command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecConcurrentOpt);

/** dialog execution method for the presolve command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecPresolve);

/** dialog execution method for the quit command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecQuit);

/** dialog execution method for the read command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecRead);

/** dialog execution method for the set default command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDefault);

/** dialog execution method for the set load command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLoad);

/** dialog execution method for the set save command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSave);

/** dialog execution method for the set diffsave command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDiffsave);

/** dialog execution method for the set parameter command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetParam);

/** dialog description method for the set parameter command */
SCIP_EXPORT
SCIP_DECL_DIALOGDESC(SCIPdialogDescSetParam);

/** dialog execution method for the fix parameter command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecFixParam);

/** dialog description method for the fix parameter command */
SCIP_EXPORT
SCIP_DECL_DIALOGDESC(SCIPdialogDescFixParam);

/** dialog execution method for the set branching direction command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingDirection);

/** dialog execution method for the set branching priority command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingPriority);

/** dialog execution method for the set heuristics aggressive command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsAggressive);

/** dialog execution method for the set heuristics default command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsDefault);

/** dialog execution method for the set heuristics fast command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsFast);

/** dialog execution method for the set heuristics off command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsOff);

/** dialog execution method for the set presolving aggressive command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingAggressive);

/** dialog execution method for the set presolving default command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingDefault);

/** dialog execution method for the set presolving fast command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingFast);

/** dialog execution method for the set presolving off command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingOff);

/** dialog execution method for the set separating aggressive command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingAggressive);

/** dialog execution method for the set separating default command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingDefault);

/** dialog execution method for the set separating fast command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingFast);

/** dialog execution method for the set separating off command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingOff);

/** dialog execution method for the set emphasis counter command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCounter);

/** dialog execution method for the set emphasis cpsolver command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCpsolver);

/** dialog execution method for the set emphasis easy CIP command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisEasycip);

/** dialog execution method for the set emphasis feasibility command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisFeasibility);

/** dialog execution method for the set emphasis hard LP command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisHardlp);

/** dialog execution method for the set emphasis optimality command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisOptimality);

/** dialog execution method for the set emphasis numerics command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisNumerics);

/** dialog execution method for the set emphasis benchmark command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisBenchmark);

/** dialog execution method for the set limits objective command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLimitsObjective);

/** dialog execution method for linear constraint type classification */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayLinearConsClassification);

/** creates a root dialog */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   );

/** @} */

/**@addtogroup DialogIncludes
 *
 * @{
 */

/** includes or updates the default dialog menus in SCIP except for menus "fix" and "set" */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeDialogDefaultBasic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes or updates the "set" menu for each available parameter setting */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeDialogDefaultSet(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes or updates the "fix" menu for each available parameter setting */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeDialogDefaultFix(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
