#ifndef _BASIC_CANVAS_H_
#define _BASIC_CANVAS_H_

#include "basic_common.h"
#include <string>

class Map;
class PointSet;

class CameraBase {
public:
	virtual float pixelGLRatio(float x, float y, float z) const = 0;
};

class Canvas 
{
public:
	Canvas() {}
	virtual ~Canvas() {}

	virtual void update_graphics() = 0;
	virtual void update_all() = 0;

	virtual CameraBase* get_camera() const = 0;
};


#endif