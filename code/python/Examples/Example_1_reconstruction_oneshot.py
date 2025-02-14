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


"""
This example demonstrates how to perform polygonal surface reconstruction
using a single function call. This approach is straightforward to use, but
running the reconstruction multiple times with different parameters will
take longer since intermediate results are not accessible.

For a step-by-step reconstruction where intermediate results can be accessed,
refer to "Example_2".
"""

import sys
sys.path.append("../../../cmake-build-release/lib/python")     # <--- Update this to use your actual build directory.

import polyfit


# Initialize polyfit
polyfit.initialize()

input_file  = "../../../data/toy_data.bvg"              # <--- Update this path to your point cloud
output_file = "../../../data/toy_data-result.obj"       # <--- Update this path to where your want to save the result

# Load the point cloud from the input file
point_cloud = polyfit.read_point_set(input_file)
if not point_cloud:
    print(f"Failed loading point cloud from file: {input_file}", file=sys.stderr)
    sys.exit(1)

# Perform the reconstruction
mesh = polyfit.reconstruct(
    point_cloud,  # input point cloud
    polyfit.SCIP, # solver name (GUROBI requires a license, which is free for research and eduction)
    0.43,  # Weight of data fitting           <--- tune if needed.
    0.27,  # Weight of model coverage         <--- tune if needed.
    0.3,   # Weight of model complexity       <--- tune if needed.
)
if not mesh:
    print("Reconstruction failed", file=sys.stderr)
    sys.exit(1)

# Output the number of faces in the reconstructed mesh
print(f"Reconstructed mesh has {mesh.size_of_facets()} faces")
# Save the reconstructed mesh to the output file
if polyfit.save_mesh(output_file, mesh):
    print(f"Reconstructed model saved to file: {output_file}")
else:
    print(f"Failed saving reconstructed model to file: {output_file}", file=sys.stderr)