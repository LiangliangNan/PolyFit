#ifndef POLYGONAL_SURFACE_RECONSTRUCTION_H
#define POLYGONAL_SURFACE_RECONSTRUCTION_H

#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>


namespace easy3d {

    namespace internal {
        class Hypothesis;
    }

    /*!
    \brief

    Implementation of PolyFit~\cite{nan2017polyfit}.
    Given a point cloud of a piecewise planar object and the initial set of planar primitives,
    PolyFit outputs a simplified and watertight surface mesh interpolating the input point cloud.

    The method first generates a set of face candidates by intersecting the planar
    primitives. Then an optimal subset of the candidate faces is selected through
    optimization under hard constraints that enforce the final model to be manifold
    and watertight.

    The reconstruction assumes the planar segmentation of the point cloud is provided in the input.
    */
    class PolyFit {
    public:

        /*!
        Creates a PolyFit object.
        After construction, candidate faces are generated and point/face confidence values are
        computed, allowing to reuse them in the subsequent reconstruction step with different parameters.
        */
        PolyFit(const PointCloud* cloud, PointCloud::VertexProperty<int> plane_indices);
        ~PolyFit();


        /** Reconstructs a watertight polygonal mesh model.
        \return the reconstructed surface mesh (empty if the reconstruction failed).
        */
        SurfaceMesh* reconstruct(
                double wt_fitting = 0.43,     ///< weight for the data fitting term.
                double wt_coverage = 0.27,    ///< weight for the point coverage term.
                double wt_complexity = 0.30   ///< weight for the model complexity term.
        );

        /*! Gives the user the possibility to access the intermediate candidate faces
                (i.e., the faces induced by the intersection of the supporting planes).
        */
        const SurfaceMesh* candidate_faces() const { return candidate_faces_; }

        /// Gets the error message (if reconstruction failed).
        const std::string &error_message() const { return error_message_; }

        // Data members.
    private:
        internal::Hypothesis* hypothesis_;

        // The generated candidate faces stored as a polygon mesh
        SurfaceMesh* candidate_faces_;

        std::string error_message_;

    private: // Copying is not allowed
        PolyFit(const PolyFit &psr);

    };

} // namespace easy3d

#endif // POLYGONAL_SURFACE_RECONSTRUCTION_H
