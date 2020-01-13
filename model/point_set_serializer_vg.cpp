/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


#include "point_set_serializer_vg.h"
#include "../basic/basic_types.h"
#include "../basic/logger.h"
#include "../basic/progress.h"
#include "../basic/color.h"
#include "../model/point_set.h"
#include "../model/vertex_group.h"
#include "../model/iterators.h"
#include <cassert>


/*
// file format definition
num_points: num
x  y  z
...

num_colors: num
r g b
...

num_normals: num
nx  ny  nz

num_groups: num

group_type: type (integer: 	VG_PLANE = 0, VG_CYLINDER = 1, VG_SPHERE = 2, VG_CONE = 3, VG_TORUS = 4, VG_GENERAL = 5)
num_group_parameters: NUM_GROUP_PARAMETERS   // number of floating point values (integer)
group_parameters: float[NUM_GROUP_PARAMETERS]
group_label: label  // the first group info
group_color: color (r, g, b)
group_num_points: num
idx ...

num_children: num

group_type: type (integer: 	VG_PLANE = 0, VG_CYLINDER = 1, VG_SPHERE = 2, VG_CONE = 3, VG_TORUS = 4, VG_GENERAL = 5)
num_group_parameters: NUM_GROUP_PARAMETERS   // number of floating point values (integer)
group_parameters: float[NUM_GROUP_PARAMETERS]
group_label: label  // 0th child of group 0
group_color: color (r, g, b)
group_num_points: num
idx ...

group_type: type (integer: 	VG_PLANE = 0, VG_CYLINDER = 1, VG_SPHERE = 2, VG_CONE = 3, VG_TORUS = 4, VG_GENERAL = 5)
num_group_parameters: NUM_GROUP_PARAMETERS   // number of floating point values (integer)
group_parameters: float[NUM_GROUP_PARAMETERS]
group_label: label  // 1st child of group 0
group_color: color (r, g, b)
group_num_points: num
idx ...
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
group_type: type (integer: 	VG_PLANE = 0, VG_CYLINDER = 1, VG_SPHERE = 2, VG_CONE = 3, VG_TORUS = 4, VG_GENERAL = 5)
num_group_parameters: NUM_GROUP_PARAMETERS   // number of floating point values (integer)
group_parameters: float[NUM_GROUP_PARAMETERS]
group_label: label  // the first group info
group_color: color (r, g, b)
group_num_points: num
idx ...
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
	for (int i = 0; i < num; ++i) {
		input >> points[i];

		std::streamoff pos = input.tellg();
		progress.notify(pos);
	}

	input >> dumy >> num;
	std::vector<vec3>& colors = pset->colors();
	colors.resize(num);
	for (int i = 0; i < num; ++i) {
		input >> colors[i];

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
