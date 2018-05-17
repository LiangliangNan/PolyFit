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