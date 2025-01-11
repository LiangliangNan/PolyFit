/* ---------------------------------------------------------------------------
 * Copyright (C) 2017 Liangliang Nan <liangliang.nan@gmail.com>
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of PolyFit. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 *
 *     Liangliang Nan and Peter Wonka.
 *     PolyFit: Polygonal Surface Reconstruction from Point Clouds.
 *     ICCV 2017.
 *
 *  For more information:
 *  https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html
 * ---------------------------------------------------------------------------
 */


#ifndef _POLYFIT_RECONSTRUCTION_H_
#define _POLYFIT_RECONSTRUCTION_H_

#include "method_common.h"
#include "../math/linear_program_solver.h"

class Map;
class PointSet;

/**
 * @brief Reconstruct a general polygonal mesh from a point cloud using PolyFit.
 * @details This function achieves reconstruction in a single function call. If you want to access the intermediate
 *      steps, you should use the HypothesisGenerator and FaceSelection classes.
 *      Check out the examples for these two use cases.
 *      For ore technical details, please refer to the paper.
 *              Liangliang Nan and Peter Wonka.
 *              PolyFit: Polygonal Surface Reconstruction from Point Clouds.
 *              ICCV 2017.https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html
 * @param point_cloud input point cloud
 * @param solver solver name. Currently, only the Gurobi (requires a license) and SCIP solvers are provided.
 * @param data_fitting weight for data fitting term
 * @param model_coverage weight for model coverage term
 * @param model_complexity weight for model complexity term
 */
Map* METHOD_API reconstruct(
        PointSet *point_cloud,                  // input point cloud
        LinearProgramSolver::SolverName solver, // solver name
        float data_fitting = 0.43,              // weight for data fitting term
        float model_coverage = 0.27,            // weight for model coverage term
        float model_complexity = 0.3            // weight for model complexity term
);


#endif
