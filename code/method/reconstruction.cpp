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


#include <method/reconstruction.h>
#include <method/face_selection.h>
#include <method/hypothesis_generator.h>
#include <method/method_global.h>

#include <basic/logger.h>
#include <model/point_set.h>


Map* reconstruct(
        PointSet *point_cloud,                  // input point cloud
        LinearProgramSolver::SolverName solver, // solver name
        float data_fitting,                     // weight for data fitting term
        float model_coverage,                   // weight for model coverage term
        float model_complexity                  // weight for model complexity term
)
{
    // step 1: refine planes
    const std::vector<VertexGroup::Ptr>& groups = point_cloud->groups();
    if (groups.empty()) {
        std::cerr << "planar segments do not exist" << std::endl;
        return nullptr;
    }
    HypothesisGenerator hypothesis(point_cloud);
    hypothesis.refine_planes();

    // step 2: generate face hypothesis
    Map* mesh = hypothesis.generate();
    if (!mesh) {
        std::cerr << "failed generating candidate faces. Please check if the input point cloud has good planar segments" << std::endl;
        return nullptr;
    }
    hypothesis.compute_confidences(mesh, false);

    // step 3: face selection
    FaceSelection selector(point_cloud, mesh);
    selector.optimize(&hypothesis, solver, data_fitting, model_coverage, model_complexity);

    if (mesh->size_of_facets() == 0) {
        std::cerr << "optimization failed: result has no face" << std::endl;
        return nullptr;
    }
    // now we don't need the point cloud anymore, and it can be deleted
    return mesh;
}

