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


#ifndef _POINT_SERIALIZER_VERTEX_GROUP_H_
#define _POINT_SERIALIZER_VERTEX_GROUP_H_

#include "model_common.h"


#include <string>
#include <vector>


class PointSet;
class VertexGroup;

class MODEL_API PointSetSerializer_vg
{
public:
	// labeled vertex groups
	static void load_vg(PointSet* pset, const std::string& file_name);
	static void save_vg(const PointSet* pset, const std::string& file_name);

	static void load_bvg(PointSet* pset, const std::string& file_name);
	static void save_bvg(const PointSet* pset, const std::string& file_name);

private:
	static VertexGroup* read_ascii_group(std::istream& input);
	static void write_ascii_group(std::ostream& output, VertexGroup* g);

	static VertexGroup* read_binary_group(std::istream& input);
	static void write_binary_group(std::ostream& output, VertexGroup* g);

	// string are stored as array of chars in binary file
	static std::string read_binary_string(std::istream& input);
	static void write_binary_string(std::ostream& output, const std::string& str);

	static std::vector<float> get_group_parameters(VertexGroup* g);
	static void assign_group_parameters(VertexGroup* g, std::vector<float>& para);
};

#endif

