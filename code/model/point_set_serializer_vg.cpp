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


#include "point_set_serializer_vg.h"

#include <cassert>

#include "../basic/basic_types.h"
#include "../basic/logger.h"
#include "../basic/progress.h"
#include "../basic/color.h"
#include "../model/point_set.h"


//#define TRANSLATE_RELATIVE_TO_FIRST_POINT


/*
// file format definition
num_points: num
x  y  z
...

num_colors: num		// can be 0; if not, it must equal to num_points
r g b               // all values are floating point numbers in range [0,1]
...

num_normals: num	// can be 0; if not, it must equal to num_points
nx  ny  nz          // all values are floating point numbers in range [0,1]

num_groups: num     // can be 0

group_type: type    // integer: PLANE = 0, CYLINDER = 1, SPHERE = 2, CONE = 3, TORUS = 4, GENERAL = 5
num_group_parameters: NUM_GROUP_PARAMETERS      // integer: the number of parameters for this group type.
group_parameters: float[NUM_GROUP_PARAMETERS]   // NUM_GROUP_PARAMETERS values in range [0,1]. Parameters of the group.
group_label: label                              // string: the label of the first group.
group_color: color (r, g, b)                    // vec3: three values in range [0,1]. The color of the group.
group_num_points: num	                        // integer: can be 0. The number of points in the group.
idx ...                                         // integers: indices of the points in the group. Total number is group_num_points.

num_children: num   // can be 0

group_type: type
num_group_parameters: NUM_GROUP_PARAMETERS
group_parameters: float[NUM_GROUP_PARAMETERS]
group_label: label
group_color: color (r, g, b)
group_num_points: num
idx ...

group_type: type
num_group_parameters: NUM_GROUP_PARAMETERS
group_parameters: float[NUM_GROUP_PARAMETERS]
group_label: label
group_color: color (r, g, b)
group_num_points: num
idx ...

...

*/
void PointSetSerializer_vg::save_vg(const PointSet* pset, const std::string& file_name) {
	// open file
	std::ofstream output(file_name.c_str());
	if (output.fail()) {
		Logger::err("-") << "could not open file\'" << file_name << "\'" << std::endl;
		return;
	}
	output.precision(16);

	//////////////////////////////////////////////////////////////////////////

	const std::vector<vec3>& points = pset->points();
	const std::vector<vec3>& colors = pset->colors();
	const std::vector<vec3>& normals = pset->normals();
	const std::vector<VertexGroup::Ptr>& groups = pset->groups();
	ProgressLogger progress(points.size() + colors.size() + normals.size() + groups.size());

	output << "num_points: " << points.size() << std::endl;
	for (std::size_t i = 0; i < points.size(); ++i) {
		output << points[i] << " ";
		progress.next();
	}
	output << std::endl;

	output << "num_colors: " << colors.size() << std::endl;
	for (std::size_t i = 0; i < colors.size(); ++i) {
		output << colors[i] << " ";
		progress.next();
	}
	output << std::endl;

	output << "num_normals: " << normals.size() << std::endl;
	for (std::size_t i = 0; i < normals.size(); ++i) {
		output << normals[i] << " ";
		progress.next();
	}
	output << std::endl;

	output << "num_groups: " << groups.size() << std::endl;
	for (std::size_t i = 0; i < groups.size(); ++i) {
		VertexGroup* g = groups[i];
		write_ascii_group(output, g);

		// children
		const std::vector<VertexGroup*>& children = g->children();
		output << "num_children: " << children.size() << std::endl;
		for (unsigned int j = 0; j < children.size(); ++j) {
			VertexGroup* chld = children[j];
			write_ascii_group(output, chld);
		}
		progress.next();
	}
}

/*
group_type: type    // integer: PLANE = 0, CYLINDER = 1, SPHERE = 2, CONE = 3, TORUS = 4, GENERAL = 5
num_group_parameters: NUM_GROUP_PARAMETERS      // integer: the number of parameters for this group type.
group_parameters: float[NUM_GROUP_PARAMETERS]   // NUM_GROUP_PARAMETERS values in range [0,1]. Parameters of this group.
group_label: label                              // string: the label of this group.
group_color: color (r, g, b)                    // vec3: three values in range [0,1]. The color of this group.
group_num_points: num	                        // integer: can be 0. The number of points in this group.
idx ...                                         // integers: indices of the points in this group. Total number is group_num_points.
*/
void PointSetSerializer_vg::write_ascii_group(std::ostream& output, VertexGroup* g) {
	//int type = g->type();
	int type = 0;
	output << "group_type: " << type << std::endl;

	const std::vector<float>& para = get_group_parameters(g);
	output << "num_group_parameters: " << para.size() << std::endl;
	output << "group_parameters: ";
	for (std::size_t i = 0; i < para.size(); ++i)
		output << para[i] << " ";
	output << std::endl;

	std::string label = g->label();
	output << "group_label: " << label << std::endl;

	Color c = g->color();
	output << "group_color: " << c.r() << " " << c.g() << " " << c.b() << std::endl;

	std::size_t num_point = g->size();
	output << "group_num_point: " << num_point << std::endl;

	for (std::size_t i = 0; i < g->size(); ++i) {
		output << g->at(i) << " ";
	}
	output << std::endl;
}


void PointSetSerializer_vg::load_vg(PointSet* pset, const std::string& file_name) {
	std::ifstream input(file_name.c_str());
	if (input.fail()) {
		Logger::err("-") << "could not open file\'" << file_name << "\'" << std::endl;
		return;
	}

	// get length of file
	input.seekg(0, input.end);
	std::streamoff length = input.tellg();
	input.seekg(0, input.beg);
	ProgressLogger progress(length);

	std::string dumy;
	std::size_t num;

	input >> dumy >> num;
	std::vector<vec3>& points = pset->points();
	points.resize(num);

#ifdef TRANSLATE_RELATIVE_TO_FIRST_POINT
	double x0, y0, z0;
	input >> x0 >> y0 >> z0;
    points[0] = vec3(0, 0, 0);

    double x, y, z;
    for (int i = 1; i < num; ++i) {
		input >> x >> y >> z;
		points[i] = vec3(x-x0, y-y0, z-z0);
		std::streamoff pos = input.tellg();
		progress.notify(pos);
	}
#else
	for (int i = 0; i < num; ++i) {
		input >> points[i];

		std::streamoff pos = input.tellg();
		progress.notify(pos);
	}
#endif
	input >> dumy >> num;
	std::vector<vec3>& colors = pset->colors();
	colors.resize(num);
	for (int i = 0; i < num; ++i) {
		input >> colors[i];

        // in case the color values are in [0, 255], converting to [0, 1]
        if (colors[i].x > 1.0f || colors[i].y > 1.0f || colors[i].z > 1.0f)
            colors[i] /= 255.0f;

		std::streamoff pos = input.tellg();
		progress.notify(pos);
	}

	input >> dumy >> num;
	std::vector<vec3>& normals = pset->normals();
	normals.resize(num);
	for (int i = 0; i < num; ++i) {
		input >> normals[i];

		std::streamoff pos = input.tellg();
		progress.notify(pos);
	}

	//////////////////////////////////////////////////////////////////////////

	std::size_t num_groups = 0;
	input >> dumy >> num_groups;
	for (int i = 0; i<num_groups; ++i) {
		VertexGroup::Ptr g = read_ascii_group(input);
        if (!g)
            continue;

		if (!g->empty()) {
			g->set_point_set(pset);
			pset->groups().push_back(g);
		}

		int num_children = 0;
		input >> dumy >> num_children;
		for (int j = 0; j<num_children; ++j) {
			VertexGroup::Ptr chld = read_ascii_group(input);
			if (!chld->empty()) {
				chld->set_point_set(pset);
				g->add_child(chld);
			}
		}

		std::streamoff pos = input.tellg();
		progress.notify(pos);
	}
}


VertexGroup* PointSetSerializer_vg::read_ascii_group(std::istream& input) {
	std::string dumy;
	int type;
	input >> dumy >> type;

	int num;
	input >> dumy >> num;
    if (num != 4)
        return nullptr;     // bad/unknown data

	std::vector<float> para(num);
	input >> dumy;
	for (int i = 0; i < num; ++i)
		input >> para[i];

	std::string label;
	input >> dumy >> label;

	float r, g, b;
	input >> dumy >> r >> g >> b;
    // in case the color values are in [0, 255], converting to [0, 1]
    if (r > 1.0f || g > 1.0f || b > 1.0f) {
        r /= 255.0f; g /= 255.0f; b /= 255.0f;
    }
	Color color(r, g, b);

	int num_points;
	input >> dumy >> num_points;

	VertexGroup* grp = new VertexGroup;
	assign_group_parameters(grp, para);

	for (int i = 0; i < num_points; ++i) {
		int idx;
		input >> idx;
		grp->push_back(idx);
	}

	grp->set_label(label);
	grp->set_color(color);

	return grp;
}


void PointSetSerializer_vg::load_bvg(PointSet* pset, const std::string& file_name) {
	std::ifstream input(file_name.c_str(), std::fstream::binary);
	if (input.fail()) {
		Logger::err("-") << "could not open file\'" << file_name << "\'" << std::endl;
		return;
	}

	int num;
	input.read((char*)(&num), sizeof(int));
	if (num <= 0) {
		Logger::err("-") << "no point exists in file\'" << file_name << "\'" << std::endl;
		return;
	}
	// read the points block
	std::vector<vec3>& points = pset->points();
	points.resize(num);
	input.read((char*)points.data(), num * sizeof(vec3));

	// read the colors block if exists
	input.read((char*)(&num), sizeof(int));
	if (num > 0) {
		if (num != points.size()) {
			Logger::err("-") << "color-point number not match" << std::endl;
			return;
		}
		std::vector<vec3>& colors = pset->colors();
		colors.resize(num);
		input.read((char*)colors.data(), num * sizeof(vec3));
	}

	// read the normals block if exists
	input.read((char*)(&num), sizeof(int));
	if (num > 0) {
		if (num != points.size()) {
			Logger::err("-") << "normal-point number not match" << std::endl;
			return;
		}
		std::vector<vec3>& normals = pset->normals();
		normals.resize(num);
		input.read((char*)normals.data(), num * sizeof(vec3));
	}

	//////////////////////////////////////////////////////////////////////////

	int num_groups = 0;
	input.read((char*)&num_groups, sizeof(int));
	for (int i = 0; i < num_groups; ++i) {
		VertexGroup::Ptr g = read_binary_group(input);
        if (!g)
            continue;

		if (!g->empty()) {
			g->set_point_set(pset);
			pset->groups().push_back(g);
		}

		int num_children = 0;
		input.read((char*)&num_children, sizeof(int));
		for (int j = 0; j < num_children; ++j) {
			VertexGroup::Ptr chld = read_binary_group(input);
			if (!chld->empty()) {
				chld->set_point_set(pset);
				g->add_child(chld);
			}
			else
				delete chld;
		}
	}
}


void PointSetSerializer_vg::save_bvg(const PointSet* pset, const std::string& file_name) {
	// open file
	std::ofstream output(file_name.c_str(), std::fstream::binary);
	if (output.fail()) {
		Logger::err("-") << "could not open file\'" << file_name << "\'" << std::endl;
		return;
	}

	// write the points block
	const std::vector<vec3>& points = pset->points();
	std::size_t num = points.size();
	output.write((char*)&num, sizeof(int));
	output.write((char*)points.data(), num * sizeof(vec3));

	const std::vector<vec3>& colors = pset->colors();
	num = colors.size();
	output.write((char*)&num, sizeof(int));
	if (num > 0)
		output.write((char*)colors.data(), num * sizeof(vec3));

	const std::vector<vec3>& normals = pset->normals();
	num = normals.size();
	output.write((char*)&num, sizeof(int));
	if (num > 0)
		output.write((char*)normals.data(), num * sizeof(vec3));

	//////////////////////////////////////////////////////////////////////////

	const std::vector<VertexGroup::Ptr>& groups = pset->groups();
	std::size_t num_groups = groups.size();
	output.write((char*)&num_groups, sizeof(int));

	for (std::size_t i = 0; i < num_groups; ++i) {
		VertexGroup* g = groups[i];
		write_binary_group(output, g);

		// children
		const std::vector<VertexGroup*>& children = g->children();
		std::size_t chld_num = children.size();
		output.write((char*)&chld_num, sizeof(int));
		for (std::size_t j = 0; j < chld_num; ++j) {
			VertexGroup* chld = children[j];
			write_binary_group(output, chld);
		}
	}
}


// string are stored as array of chars in binary file
std::string PointSetSerializer_vg::read_binary_string(std::istream& input) {
	int size_int = sizeof(int);
	int size_char = sizeof(char);

	int num_char = 0;
	input.read((char*)&num_char, size_int);

	std::string str;
	for (int i = 0; i < num_char; ++i) {
		char c;
		input.read(&c, size_char);
		if (c != ' ')
			str.push_back(c);
		else
			str.push_back('-');
	}

	return str;
}


// string are stored as array of chars in binary file
void PointSetSerializer_vg::write_binary_string(std::ostream& output, const std::string& str) {
	int size_int = sizeof(int);
	int size_char = sizeof(char);

	std::size_t num_char = str.size();
	output.write((char*)&num_char, size_int);

	for (std::size_t i = 0; i < num_char; ++i) {
		char c = str[i];
		if (c == ' ')
			c = '-';
		output.write(&c, size_char);
	}
}


// for binary file format, no string stuff except labels. we add size info before each label
VertexGroup* PointSetSerializer_vg::read_binary_group(std::istream& input) {
	int type;
	input.read((char*)&type, sizeof(int));

	int num;
	input.read((char*)&num, sizeof(int));
    if (num != 4)
        return nullptr;     // bad/unknown data

	std::vector<float> para(num);
	input.read((char*)para.data(), num * sizeof(float));

	VertexGroup* grp = new VertexGroup;
	assign_group_parameters(grp, para);

	////////////////////////////////////////////////////////////////////////// 
	// Liangliang: check why this doesn't work in release mode
	//int label_size = 0;   
	//std::string label;
	//input.read((char*)&label_size, sizeof(int));  
	//input.read((char*)&label, label_size); 
	//////////////////////////////////////////////////////////////////////////
	std::string label = read_binary_string(input);
	grp->set_label(label);

	float arr[3];
	input.read((char*)arr, 3 * sizeof(float));
	Color color(arr);
	grp->set_color(color);

	int num_points = 0;
	input.read((char*)&num_points, sizeof(int));
	grp->resize(num_points);
	input.read((char*)grp->data(), num_points * sizeof(int));

	return grp;
}


// for binary file format, no string stuff except labels. we add size info before each label
void PointSetSerializer_vg::write_binary_group(std::ostream& output, VertexGroup* g) {
	//int type = g->type();
	int type = 0;
	output.write((char*)&type, sizeof(int));

	int num = 4;
	output.write((char*)&num, sizeof(int));

	const std::vector<float>& para = get_group_parameters(g);
	output.write((char*)para.data(), sizeof(float) * num);

	//////////////////////////////////////////////////////////////////////////
	// Liangliang: check why not work
	// 	std::string label = g->label();
	// 	int label_size = sizeof(label) + 1;
	// 	output.write((char*)&label_size, sizeof(int));
	// 	output.write(label.c_str(), label_size + 1);
	//////////////////////////////////////////////////////////////////////////
	std::string label = g->label();
	write_binary_string(output, label);

	Color c = g->color();
	output.write((char*)c.data(), sizeof(float) * 3);

	std::size_t num_point = g->size();
	output.write((char*)&num_point, sizeof(int));
	output.write((char*)g->data(), num_point * sizeof(int));
}


std::vector<float> PointSetSerializer_vg::get_group_parameters(VertexGroup* g) {
	int num = 4;
	std::vector<float> para(num);

	para[0] = static_cast<float>(g->plane().a());
	para[1] = static_cast<float>(g->plane().b());
	para[2] = static_cast<float>(g->plane().c());
	para[3] = static_cast<float>(g->plane().d());

	return para;
}


void PointSetSerializer_vg::assign_group_parameters(VertexGroup* g, std::vector<float>& para) {
	int num = 4;
	assert(para.size() == num);

	g->set_plane(Plane3d(para[0], para[1], para[2], para[3]));
}
