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

