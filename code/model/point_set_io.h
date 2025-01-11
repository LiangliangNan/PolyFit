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

#ifndef _POINT_SET_IO_H_
#define _POINT_SET_IO_H_

#include "model_common.h"

#include <string>


class PointSet;

class MODEL_API PointSetIO
{
public:
	// for both point cloud and mesh
	static PointSet* read(const std::string& file_name);

	// save the point set to a file. return false if failed.
	static bool		 save(const std::string& file_name, const PointSet* point_set);
};

#endif