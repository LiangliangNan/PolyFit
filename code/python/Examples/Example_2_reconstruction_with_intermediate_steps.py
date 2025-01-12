# ------------------------------------------------------------------
# Copyright (C) 2025 Liangliang Nan <liangliang.nan@gmail.com>
# https://3d.bk.tudelft.nl/liangliang/
#
# This file is part of PolyFit. If it is useful in your research/work,
# I would be grateful if you show your appreciation by citing it:
#
#     Liangliang Nan and Peter Wonka.
#     PolyFit: Polygonal Surface Reconstruction from Point Clouds.
#     ICCV 2017.
# ------------------------------------------------------------------


# This example shows how to perform the reconstruction step by step.
# This could be useful if you want to reuse the generated candidate faces and their confidence values,
# but try multiple times of the final optimization step with different parameters.
#
# See "Example_1" for reconstruction with a single function call.


import sys
sys.path.append("../../../cmake-build-release/lib/python")     # <--- Update this to use your actual build directory.

import PyPolyFit as polyfit

def main():
    # Initialize PolyFit
    polyfit.initialize()

    input_file  = "../../../data/toy_data.bvg"          # <--- Update this path to your point cloud
    output_file = "../../../data/toy_data-result.obj"   # <--- Update this path to where your want to save the result

    # Load point cloud from file
    point_cloud = polyfit.read_point_set(input_file)
    if not point_cloud:
        print(f"Failed loading point cloud from file: {input_file}", file=sys.stderr)
        sys.exit(1)

    # Step 1: Refine planes
    print("Refining planes...")
    groups = point_cloud.groups()
    if not groups:
        print("Planar segments do not exist. Read the ReadMe file carefully!!!", file=sys.stderr)
        sys.exit(1)

    # Create a hypothesis generator
    hypothesis = polyfit.HypothesisGenerator(point_cloud)
    hypothesis.refine_planes()      # refine planes

    # Step 2: Generate face hypothesis
    print("Generating face hypothesis...")
    mesh = hypothesis.generate()    # generate candidate faces
    if not mesh:
        print("Failed generating candidate faces. Please check if the input point cloud has good planar segments", file=sys.stderr)
        sys.exit(1)

    # Compute confidences of candidate faces
    hypothesis.compute_confidences(mesh, False)

    # Step 3: Face selection
    print("Optimization...")
    # Create an instance of FaceSelection
    selector = polyfit.FaceSelection(point_cloud, mesh)
    selector.optimize(hypothesis,
                      polyfit.SCIP, # solver name (GUROBI requires a license, free for research and eduction)
                      0.43,         # Weight of data_fitting         <--- tune if needed.
                      0.27,         # Weight of model_coverage       <--- tune if needed.
                      0.3)          # Weight of model_complexity     <--- tune if needed.

    if mesh.size_of_facets() == 0:
        print("Optimization failed: model has no face", file=sys.stderr)
        sys.exit(1)
    print(f"Reconstructed mesh has {mesh.size_of_facets()} faces")

    # Save result to file
    if polyfit.save_mesh(output_file, mesh):
        print(f"Reconstructed model saved to file: {output_file}")
    else:
        print(f"Failed saving reconstructed model to file: {output_file}", file=sys.stderr)
        sys.exit(1)

    sys.exit(0)

if __name__ == "__main__":
    main()
