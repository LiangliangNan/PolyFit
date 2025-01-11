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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "../../basic/logger.h"
#include "../../model/point_set.h"
#include "../../model/map.h"
#include "../../method/method_global.h"
#include "../../method/hypothesis_generator.h"
#include "../../method/face_selection.h"
#include "../../method/reconstruction.h"
#include "../../model/map_io.h"
#include "../../model/point_set_io.h"
#include "../../math/linear_program_solver.h"


namespace py = pybind11;

// Define Python bindings for the SmartPointer template class
template <typename T>
void bind_smart_pointer(py::module &m, const std::string &name) {
    py::class_<SmartPointer<T>>(m, name.c_str())
            .def(py::init<>())
            .def(py::init<T*>())
            .def(py::init<const SmartPointer<T>&>())
            .def("get", (T* (SmartPointer<T>::*)()) &SmartPointer<T>::get, py::return_value_policy::reference)
            .def("forget", &SmartPointer<T>::forget)
            .def("is_nil", &SmartPointer<T>::is_nil)
            .def("__eq__", [](const SmartPointer<T>& self, T* ptr) { return self.get() == ptr; })
            .def("__ne__", [](const SmartPointer<T>& self, T* ptr) { return self.get() != ptr; })
            .def("__bool__", [](const SmartPointer<T>& self) { return !self.is_nil(); })
            .def("__repr__", [](const SmartPointer<T>& self) {
                return "<SmartPointer to " + std::to_string(reinterpret_cast<uintptr_t>(self.get())) + ">";
            });
}

// Define Python bindings for the VertexGroup class
void bind_vertex_group(py::module &m) {
    // Bind the SmartPointer<VertexGroup> type
    bind_smart_pointer<VertexGroup>(m, "SmartPointerVertexGroup");

    // Bind the VertexGroup class
    py::class_<VertexGroup>(m, "VertexGroup")
            .def("nb_vertice", &VertexGroup::nb_vertice, "Get the number of vertices in the group")
            .def("set_point_set", &VertexGroup::set_point_set, "Set the associated point set")
            .def("plane", &VertexGroup::plane, "Get the plane of the vertex group")
            .def("set_plane", &VertexGroup::set_plane, "Set the plane of the vertex group");
}


// Define Python bindings for the PointSet class
void bind_point_set(py::module &m) {
    py::class_<PointSet>(m, "PointSet")
            .def(py::init<>())
            .def("num_points", &PointSet::num_points, "Get the number of points")
            .def("points", py::overload_cast<>(&PointSet::points), "Get the list of points")
            .def("colors", py::overload_cast<>(&PointSet::colors), "Get the list of colors")
            .def("normals", py::overload_cast<>(&PointSet::normals), "Get the list of normals")
            .def("planar_qualities", py::overload_cast<>(&PointSet::planar_qualities),
                 "Get the list of planar qualities")
            .def("has_normals", &PointSet::has_normals, "Check if normals are available")
            .def("has_colors", &PointSet::has_colors, "Check if colors are available")
            .def("has_planar_qualities", &PointSet::has_planar_qualities, "Check if planar qualities are available")
            .def("groups", py::overload_cast<>(&PointSet::groups), "Get the list of vertex groups")
            .def("fit_plane", &PointSet::fit_plane, "Fit a plane for a vertex group");
}


// Define Python bindings for the Map class and related classes
void bind_map(py::module &m) {
    // Bind the Map class
    py::class_<Map>(m, "Map")
            .def(py::init<>())
            .def("size_of_vertices", &Map::size_of_vertices, "Get the number of vertices")
            .def("size_of_halfedges", &Map::size_of_halfedges, "Get the number of halfedges")
            .def("size_of_facets", &Map::size_of_facets, "Get the number of facets")
            .def("is_triangulated", &Map::is_triangulated, "Check if the map is triangulated");
}


// Define Python bindings for the global variables and constants
void bind_method_global(py::module &m) {
    // Expose global variables
    m.attr("lambda_data_fitting") = &Method::lambda_data_fitting;
    m.attr("lambda_model_coverage") = &Method::lambda_model_coverage;
    m.attr("lambda_model_complexity") = &Method::lambda_model_complexity;
    m.attr("snap_sqr_distance_threshold") = &Method::snap_sqr_distance_threshold;

    // Expose global constants
    m.attr("facet_attrib_supporting_vertex_group") = Method::facet_attrib_supporting_vertex_group;
    m.attr("facet_attrib_supporting_point_num") = Method::facet_attrib_supporting_point_num;
    m.attr("facet_attrib_facet_area") = Method::facet_attrib_facet_area;
    m.attr("facet_attrib_covered_area") = Method::facet_attrib_covered_area;
}


// Define Python bindings for the HypothesisGenerator class
void bind_hypothesis_generator(py::module &m) {
    // Bind the Adjacency type
    py::bind_vector<HypothesisGenerator::Adjacency>(m, "Adjacency");

    // Bind the HypothesisGenerator class
    py::class_<HypothesisGenerator>(m, "HypothesisGenerator")
            .def(py::init<PointSet *>(), py::arg("pset"))
            .def("refine_planes", &HypothesisGenerator::refine_planes, "Refine planes")
            .def("generate", &HypothesisGenerator::generate, "Generate candidate faces")
            .def("compute_confidences", &HypothesisGenerator::compute_confidences,
                 py::arg("mesh"), py::arg("use_conficence") = false, "Compute confidences")
            .def("extract_adjacency", &HypothesisGenerator::extract_adjacency, "Extract adjacency information")
            .def("ready_for_optimization", &HypothesisGenerator::ready_for_optimization,
                 "Check if ready for optimization");
}


// Define Python bindings for the FaceSelection class
void bind_face_selection(py::module &m) {
    // Bind the SolverName enum
    py::enum_<LinearProgramSolver::SolverName>(m, "")
            .value("SCIP", LinearProgramSolver::SCIP)
            .export_values();

    // Bind the FaceSelection class
    py::class_<FaceSelection>(m, "FaceSelection")
            .def(py::init<PointSet *, Map *>(), py::arg("pset"), py::arg("model"))
            .def("optimize", &FaceSelection::optimize,
                 py::arg("generator"),
                 py::arg("solver_name") = LinearProgramSolver::SCIP,
                 "Optimize face selection"
            );
}


// Define Python bindings for the reconstruct function
void bind_reconstruction(py::module &m) {
    m.def("reconstruct", &reconstruct,
          py::arg("point_cloud"),
          py::arg("solver") = LinearProgramSolver::SolverName::SCIP,
          py::arg("data_fitting") = 0.43f,
          py::arg("model_coverage") = 0.27f,
          py::arg("model_complexity") = 0.3f,
          R"pbdoc(
            Reconstruct a general polygonal mesh from a point cloud using PolyFit.

            Parameters:
                point_cloud (PointSet): Input point cloud.
                solver (LinearProgramSolver.SolverName): Solver name (e.g., SCIP or GUROBI).
                data_fitting (float): Weight for data fitting term (default: 0.43).
                model_coverage (float): Weight for model coverage term (default: 0.27).
                model_complexity (float): Weight for model complexity term (default: 0.3).

            Returns:
                Map: Reconstructed polygonal mesh.
          )pbdoc");
}

// Define Python bindings
PYBIND11_MODULE(PyPolyFit_NAME, m) {
    m.doc() = "Python bindings for PolyFit";

    m.def("initialize", &Logger::initialize, "Initialize PolyFit");

    bind_vertex_group(m);
    bind_point_set(m);
    bind_map(m);
    bind_hypothesis_generator(m);
    bind_method_global(m);
    bind_face_selection(m);
    bind_reconstruction(m);

    // Bind I/O functions
    m.def("save_mesh", &MapIO::save, "Save a mesh to a file", py::arg("mesh"), py::arg("filename"));
    m.def("read_point_set", &PointSetIO::read, "Load a point cloud from a file", py::arg("filename"));
}