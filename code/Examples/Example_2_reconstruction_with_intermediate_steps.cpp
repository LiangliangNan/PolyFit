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


#include "../basic/logger.h"
#include "../model/point_set.h"
#include "../model/map.h"
#include "../method/method_global.h"
#include "../method/hypothesis_generator.h"
#include "../method/face_selection.h"
#include "../model/map_io.h"
#include "../model/point_set_io.h"


int main(int argc, char **argv)
{
    // initialize the logger (this is not optional)
    Logger::initialize();

    // input point cloud file name
    const std::string input_file = (argc > 1) ? argv[1] : std::string(POLYFIT_CODE_DIR) + "/../data/toy_data.bvg";
    // output mesh file name
    const std::string output_file = (argc > 2) ? argv[2] : std::string(POLYFIT_CODE_DIR) + "/../data/toy_data-result.obj";

    // below are the default parameters (change these when necessary)
    Method::lambda_data_fitting = 0.43;
    Method::lambda_model_coverage = 0.27;
    Method::lambda_model_complexity = 0.3;

    // load point cloud from file
    PointSet* pset = PointSetIO::read(input_file);
    if (!pset) {
        std::cerr << "failed loading point cloud from file: " << input_file << std::endl;
        return EXIT_FAILURE;
    }

    // step 1: refine planes
    std::cout << "refining planes..." << std::endl;
    const std::vector<VertexGroup::Ptr>& groups = pset->groups();
    if (groups.empty()) {
        std::cerr << "planar segments do not exist" << std::endl;
        return EXIT_FAILURE;
    }
    HypothesisGenerator hypothesis(pset);
    hypothesis.refine_planes();

    // step 2: generate face hypothesis
    std::cout << "generating plane hypothesis..." << std::endl;
    Map* mesh = hypothesis.generate();
    if (!mesh) {
        std::cerr << "failed generating candidate faces. Please check if the input point cloud has good planar segments" << std::endl;
        return EXIT_FAILURE;
    }
    hypothesis.compute_confidences(mesh, false);

    // step 3: face selection
    std::cout << "optimization..." << std::endl;
    FaceSelection selector(pset, mesh);

#ifdef HAS_GUROBI
    selector.optimize(&hypothesis, LinearProgramSolver::GUROBI);
#else
    selector.optimize(&hypothesis, LinearProgramSolver::SCIP);
#endif

    if (mesh->size_of_facets() == 0) {
        std::cerr << "optimization failed: model has no face" << std::endl;
        return EXIT_FAILURE;
    }
    // now we don't need the point cloud anymore, and it can be deleted
    delete pset;

    // step 4: save result to file
    if (MapIO::save(output_file, mesh))
        std::cout << "reconstructed model saved to file: " << output_file << std::endl;
    else {
        std::cerr << "failed saving reconstructed model to file: " << output_file << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


