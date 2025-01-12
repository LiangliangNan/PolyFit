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


#include <basic/logger.h>
#include <model/map.h>
#include <model/map_io.h>
#include <model/point_set.h>
#include <model/point_set_io.h>
#include <method/reconstruction.h>


int main(int argc, char **argv)
{
    // initialize the logger (this is not optional)
    Logger::initialize();

    // input point cloud file name
    const std::string input_file = (argc > 1) ? argv[1] : std::string(POLYFIT_CODE_DIR) + "/../data/toy_data.bvg";
    // output mesh file name
    const std::string output_file = (argc > 2) ? argv[2] : std::string(POLYFIT_CODE_DIR) + "/../data/toy_data-result.obj";

    // load point cloud from file
    auto point_cloud = PointSetIO::read(input_file);
    if (!point_cloud) {
        std::cerr << "failed loading point cloud from file: " << input_file << std::endl;
        return EXIT_FAILURE;
    }

    // reconstruction
    auto mesh = reconstruct(
            point_cloud,
            LinearProgramSolver::SCIP, // (GUROBI requires a license, which is free for research and eduction)
            0.43,  // Weight of data_fitting
            0.27,  // Weight of model_coverage
            0.3   // Weight of model_complexity
    );

    if (!mesh) {
        std::cerr << "Reconstruction failed" << std::endl;
        return EXIT_FAILURE;
    }
    // now we don't need the point cloud anymore, and it can be deleted
    delete point_cloud;

    std::cout << "Reconstructed mesh has " << mesh->size_of_facets() << " faces" << std::endl;
    if (MapIO::save(output_file, mesh)) {
        std::cout << "Reconstructed model saved to file: " << output_file << std::endl;
        delete mesh;
        return EXIT_SUCCESS;
    }
    else {
        std::cerr << "Failed saving reconstructed model to file: " << output_file << std::endl;
        return EXIT_FAILURE;
    }
}


