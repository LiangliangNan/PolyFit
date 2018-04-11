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

/**@file   scipdefplugins.h
 * @ingroup PUBLICCOREAPI
 * @brief  default SCIP plugins
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIPDEFPLUGINS_H__
#define __SCIP_SCIPDEFPLUGINS_H__

#include "scip/scip.h"

/* include header files here, such that the user only has to include
 * scipdefplugins.h
 */
#include "scip/branch_allfullstrong.h"
#include "scip/branch_cloud.h"
#include "scip/branch_distribution.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_multaggr.h"
#include "scip/branch_nodereopt.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"
#include "scip/compr_largestrepr.h"
#include "scip/compr_weakcompr.h"
#include "scip/cons_abspower.h"
#include "scip/cons_and.h"
#include "scip/cons_bivariate.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_cardinality.h"
#include "scip/cons_conjunction.h"
#include "scip/cons_countsols.h"
#include "scip/cons_cumulative.h"
#include "scip/cons_disjunction.h"
#include "scip/cons_indicator.h"
#include "scip/cons_integral.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_linking.h"
#include "scip/cons_logicor.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_or.h"
#include "scip/cons_orbisack.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_pseudoboolean.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_setppc.h"
#include "scip/cons_soc.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_superindicator.h"
#include "scip/cons_symresack.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"
#include "scip/cons_components.h"
#include "scip/disp_default.h"
#include "scip/event_solvingphase.h"
#include "scip/event_softtimelimit.h"
#include "scip/heur_actconsdiving.h"
#include "scip/heur_bound.h"
#include "scip/heur_clique.h"
#include "scip/heur_coefdiving.h"
#include "scip/heur_completesol.h"
#include "scip/heur_crossover.h"
#include "scip/heur_dins.h"
#include "scip/heur_distributiondiving.h"
#include "scip/heur_dualval.h"
#include "scip/heur_feaspump.h"
#include "scip/heur_fixandinfer.h"
#include "scip/heur_fracdiving.h"
#include "scip/heur_gins.h"
#include "scip/heur_guideddiving.h"
#include "scip/heur_indicator.h"
#include "scip/heur_intdiving.h"
#include "scip/heur_intshifting.h"
#include "scip/heur_linesearchdiving.h"
#include "scip/heur_localbranching.h"
#include "scip/heur_locks.h"
#include "scip/heur_lpface.h"
#include "scip/heur_alns.h"
#include "scip/heur_multistart.h"
#include "scip/heur_mutation.h"
#include "scip/heur_mpec.h"
#include "scip/heur_nlpdiving.h"
#include "scip/heur_objpscostdiving.h"
#include "scip/heur_octane.h"
#include "scip/heur_ofins.h"
#include "scip/heur_oneopt.h"
#include "scip/heur_pscostdiving.h"
#include "scip/heur_proximity.h"
#include "scip/heur_randrounding.h"
#include "scip/heur_rens.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_repair.h"
#include "scip/heur_rins.h"
#include "scip/heur_rootsoldiving.h"
#include "scip/heur_rounding.h"
#include "scip/heur_shiftandpropagate.h"
#include "scip/heur_shifting.h"
#include "scip/heur_simplerounding.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trivial.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heur_trysol.h"
#include "scip/heur_twoopt.h"
#include "scip/heur_undercover.h"
#include "scip/heur_vbounds.h"
#include "scip/heur_veclendiving.h"
#include "scip/heur_zeroobj.h"
#include "scip/heur_zirounding.h"
#include "scip/nodesel_bfs.h"
#include "scip/nodesel_breadthfirst.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_uct.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/presol_boundshift.h"
#include "scip/presol_convertinttobin.h"
#include "scip/presol_domcol.h"
#include "scip/presol_implfree.h"
#include "scip/presol_dualagg.h"
#include "scip/presol_dualcomp.h"
#include "scip/presol_dualinfer.h"
#include "scip/presol_gateextraction.h"
#include "scip/presol_implics.h"
#include "scip/presol_inttobinary.h"
#include "scip/presol_redvub.h"
#include "scip/presol_symbreak.h"
#include "scip/presol_qpkktref.h"
#include "scip/presol_trivial.h"
#include "scip/presol_tworowbnd.h"
#include "scip/presol_sparsify.h"
#include "scip/presol_stuffing.h"
#include "scip/presol_symmetry.h"
#include "scip/prop_dualfix.h"
#include "scip/prop_genvbounds.h"
#include "scip/prop_nlobbt.h"
#include "scip/prop_obbt.h"
#include "scip/prop_orbitalfixing.h"
#include "scip/prop_probing.h"
#include "scip/prop_pseudoobj.h"
#include "scip/prop_redcost.h"
#include "scip/prop_rootredcost.h"
#include "scip/prop_vbounds.h"
#include "scip/reader_bnd.h"
#include "scip/reader_ccg.h"
#include "scip/reader_cip.h"
#include "scip/reader_cnf.h"
#include "scip/reader_diff.h"
#include "scip/reader_fix.h"
#include "scip/reader_fzn.h"
#include "scip/reader_gms.h"
#include "scip/reader_lp.h"
#include "scip/reader_mps.h"
#include "scip/reader_mst.h"
#include "scip/reader_opb.h"
#include "scip/reader_osil.h"
#include "scip/reader_pip.h"
#include "scip/reader_ppm.h"
#include "scip/reader_pbm.h"
#include "scip/reader_rlp.h"
#include "scip/reader_sol.h"
#include "scip/reader_wbo.h"
#include "scip/reader_zpl.h"
#include "scip/sepa_eccuts.h"
#include "scip/sepa_cgmip.h"
#include "scip/sepa_clique.h"
#include "scip/sepa_closecuts.h"
#include "scip/sepa_aggregation.h"
#include "scip/sepa_convexproj.h"
#include "scip/sepa_disjunctive.h"
#include "scip/sepa_gauge.h"
#include "scip/sepa_gomory.h"
#include "scip/sepa_impliedbounds.h"
#include "scip/sepa_intobj.h"
#include "scip/sepa_mcf.h"
#include "scip/sepa_oddcycle.h"
#include "scip/sepa_rapidlearning.h"
#include "scip/sepa_strongcg.h"
#include "scip/sepa_zerohalf.h"
#include "scip/scipshell.h"
#include "scip/table_default.h"
#include "scip/concsolver_scip.h"

#include "nlpi/nlpi_ipopt.h"
#include "nlpi/nlpi_filtersqp.h"
#include "nlpi/nlpi_worhp.h"
#include "nlpi/nlpi_all.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes default SCIP plugins into SCIP */
EXTERN
SCIP_RETCODE SCIPincludeDefaultPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
