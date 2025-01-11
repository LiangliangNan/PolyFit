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

#include "point_set_io.h"
#include "point_set_serializer_vg.h"
#include "../model/point_set.h"
#include "../basic/stop_watch.h"
#include "../basic/file_utils.h"
#include "../basic/logger.h"

#include <fstream>
#include <list>


PointSet* PointSetIO::read(const std::string& file_name)
{
	std::ifstream in(file_name.c_str()) ;
	if(in.fail()) {
		Logger::err("-") << "cannot open file: " << file_name << std::endl;
		return nil;
	}
	in.close();

	Logger::out("-") << "reading file..." << std::endl;

	std::string ext = FileUtils::extension(file_name);
	String::to_lowercase(ext);

	StopWatch w;
	PointSet* pset = new PointSet;
	if (ext == "vg")
		PointSetSerializer_vg::load_vg(pset, file_name);
	else if (ext == "bvg")
		PointSetSerializer_vg::load_bvg(pset, file_name);

	else {
		Logger::err("-") << "reading file failed (unknown file format)" << std::endl;
		delete pset;
		return nil;
	}
		
	if (pset->num_points() < 1) {
		Logger::err("-") << "reading file failed (no data exist)" << std::endl;
		delete pset;
		return nil;
	}

	Logger::out("-") << "done. " << w.elapsed() << " sec." << std::endl;

	return pset;
}

bool PointSetIO::save(const std::string& file_name, const PointSet* point_set) {
	if (!point_set) {
		Logger::err("-") << "Point set is null" << std::endl;
		return false;
	}
	
	std::ofstream out(file_name.c_str()) ;
	if(out.fail()) {
		Logger::err("-") << "cannot open file: \'" << file_name << "\' for writing" << std::endl;
		return false;
	}
	Logger::out("-") << "saving file..." << std::endl;
	out.close();

	StopWatch w;
	std::string ext = FileUtils::extension(file_name);
	String::to_lowercase(ext);
	
	out.precision(16);

 	if (ext == "vg")
		PointSetSerializer_vg::save_vg(point_set, file_name);
	else if (ext == "bvg")
		PointSetSerializer_vg::save_bvg(point_set, file_name);

	else {
		Logger::err("-") << "saving file failed (unknown file format)" << std::endl;
		return false;
	}

	Logger::out("-") << "done. " << w.elapsed() << " sec." << std::endl;
	return true;
}