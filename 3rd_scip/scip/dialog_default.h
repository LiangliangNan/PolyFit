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

/**@file   dialog_default.h
 * @ingroup DIALOGS
 * @brief  default user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DIALOG_DEFAULT_H__
#define __SCIP_DIALOG_DEFAULT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup DIALOGS
 *
 * @{
 */

/** standard menu dialog execution method, that displays it's help screen if the remaining command line is empty */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenu);

/** standard menu dialog execution method, that doesn't display it's help screen */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenuLazy);

/** dialog execution method for the change add constraint */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeAddCons);

/** dialog execution method for the change bounds command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeBounds);

/** dialog execution method for the freetransproblem command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeFreetransproblem);

/** dialog execution method for the changing the objective sense */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeObjSense);

/** dialog execution method for the checksol command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChecksol);

/** dialog execution method for the cliquegraph command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCliquegraph);

/** dialog execution method for the display branching command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayBranching);

/** dialog execution method for the display compression command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayCompression);

/** dialog execution method for the display conflict command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConflict);

/** dialog execution method for the display conshdlrs command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConshdlrs);

/** dialog execution method for the display displaycols command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayDisplaycols);

/** dialog execution method for the display heuristics command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayHeuristics);

/** dialog execution method for the display memory command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayMemory);

/** dialog execution method for the display nodeselectors command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNodeselectors);

/** dialog execution method for the display nlpi command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNlpi);

/** dialog execution method for the display parameters command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayParameters);

/** dialog execution method for the display presolvers command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPresolvers);

/** dialog execution method for the display pricer command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPricers);

/** dialog execution method for the display problem command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayProblem);

/** dialog execution method for the display propagators command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPropagators);

/** dialog execution method for the display readers command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayReaders);

/** dialog execution method for the display relaxators command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayRelaxators);

/** dialog execution method for the display separators command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySeparators);

/** dialog execution method for the display solution command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySolution);

/** dialog execution method for the display finitesolution command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayFiniteSolution);

/** dialog execution method for the display dual solution command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayDualSolution);

/** dialog execution method for the display of solutions in the pool command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySolutionPool);

/** dialog execution method for the display statistics command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayStatistics);

/** dialog execution method for the display reoptstatistics command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayReoptStatistics);

/** dialog execution method for the display transproblem command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTransproblem);

/** dialog execution method for the display value command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayValue);

/** dialog execution method for the display varbranchstatistics command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayVarbranchstatistics);

/** dialog execution method for the display LP solution quality command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayLPSolutionQuality);

/** dialog execution method for the display transsolution command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTranssolution);

/** dialog execution method for the help command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecHelp);

/** dialog execution method for the free command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecFree);

/** dialog execution method for the newstart command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecNewstart);

/** dialog execution method for the transform command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecTransform);

/** dialog execution method for the optimize command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecOptimize);

/** dialog execution method for the parallelopt command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecConcurrentOpt);

/** dialog execution method for the presolve command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecPresolve);

/** dialog execution method for the quit command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecQuit);

/** dialog execution method for the read command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecRead);

/** dialog execution method for the set default command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDefault);

/** dialog execution method for the set load command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLoad);

/** dialog execution method for the set save command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSave);

/** dialog execution method for the set diffsave command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDiffsave);

/** dialog execution method for the set parameter command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetParam);

/** dialog description method for the set parameter command */
EXTERN
SCIP_DECL_DIALOGDESC(SCIPdialogDescSetParam);

/** dialog execution method for the fix parameter command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecFixParam);

/** dialog description method for the fix parameter command */
EXTERN
SCIP_DECL_DIALOGDESC(SCIPdialogDescFixParam);

/** dialog execution method for the set branching direction command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingDirection);

/** dialog execution method for the set branching priority command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingPriority);

/** dialog execution method for the set heuristics aggressive command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsAggressive);

/** dialog execution method for the set heuristics default command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsDefault);

/** dialog execution method for the set heuristics fast command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsFast);

/** dialog execution method for the set heuristics off command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsOff);

/** dialog execution method for the set presolving aggressive command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingAggressive);

/** dialog execution method for the set presolving default command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingDefault);

/** dialog execution method for the set presolving fast command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingFast);

/** dialog execution method for the set presolving off command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingOff);

/** dialog execution method for the set separating aggressive command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingAggressive);

/** dialog execution method for the set separating default command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingDefault);

/** dialog execution method for the set separating fast command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingFast);

/** dialog execution method for the set separating off command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingOff);

/** dialog execution method for the set emphasis counter command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCounter);

/** dialog execution method for the set emphasis cpsolver command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCpsolver);

/** dialog execution method for the set emphasis easy CIP command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisEasycip);

/** dialog execution method for the set emphasis feasibility command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisFeasibility);

/** dialog execution method for the set emphasis hard LP command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisHardlp);

/** dialog execution method for the set emphasis optimality command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisOptimality);

/** dialog execution method for the set limits objective command */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLimitsObjective);

/** dialog execution method for linear constraint type classification */
EXTERN
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayLinearConsClassification);

/** creates a root dialog */
EXTERN
SCIP_RETCODE SCIPcreateRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   );

/* @} */

/**@addtogroup DialogIncludes
 *
 * @{
 */

/** includes or updates the default dialog menus in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeDialogDefault(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes or updates the "set" menu for each available parameter setting */
EXTERN
SCIP_RETCODE SCIPincludeDialogDefaultSet(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes or updates the "fix" menu for each available parameter setting */
EXTERN
SCIP_RETCODE SCIPincludeDialogDefaultFix(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
