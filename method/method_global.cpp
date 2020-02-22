#include "method_global.h"


namespace Method {

	double lambda_data_fitting = 0.43;
	double lambda_model_coverage = 0.27;
	double lambda_model_complexity = 0.30;

	double snap_sqr_distance_threshold = 1e-14;

	//________________ names for various quality measures ____________________

	std::string facet_attrib_supporting_vertex_group = "facet_supporting_vertex_group";
	std::string facet_attrib_supporting_point_num = "facet_supporting_point_num";
	std::string facet_attrib_facet_area = "facet_area";
	std::string facet_attrib_covered_area = "facet_covered_area";


}
