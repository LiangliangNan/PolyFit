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


#ifndef _VERTEX_GROUP_H_
#define _VERTEX_GROUP_H_


#include "../basic/color.h"
#include "../math/principal_axes.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"
#include <vector>
#include <set>


class PointSet;
class VertexGroup : public std::vector<unsigned int>, public Counted
{
public:
	typedef SmartPointer<VertexGroup>	Ptr;

public:
	VertexGroup(PointSet* pset = nil) 
		: label_("unknown")
		, point_set_(pset)
		, visible_(true)
		, highlighted_(false)
	{}

	~VertexGroup() {
		clear();
	}

	const std::string& label() const { return label_; }
	void set_label(const std::string& lb) { label_ = lb; } 

	size_t nb_vertice() { return size(); }

	//////////////////////////////////////////////////////////////////////////

	PointSet* point_set() { return point_set_; }
	const PointSet* point_set() const { return point_set_; }
	void set_point_set(PointSet* pset) { point_set_ = pset; }

	const Color& color() const { return color_; }
	void set_color(const Color& c) { color_ = c; }

	void set_plane(const Plane3d& plane) { plane_ = plane; }
	const Plane3d& plane() const { return plane_; }

	const std::vector<unsigned int>& boundary() const { return boundary_; }
	void set_boundary(const std::vector<unsigned int>& bd) { boundary_ = bd; }
	
	//////////////////////////////////////////////////////////////////////////

	VertexGroup* parent() { return parent_; }
	void set_parent(VertexGroup* g) { parent_ = g; }

	std::vector<VertexGroup*> children();
	void set_children(const std::vector<VertexGroup*>& chld);

	void add_child(VertexGroup* g);
	void remove_child(VertexGroup* g);

	void remove_children();
	void delete_children();

	bool is_visible() const  { return visible_; }
	void set_visible(bool b) { visible_ = b; }

	bool is_highlighted() const  { return highlighted_; }
	virtual void set_highlighted(bool b) { highlighted_ = b; }

private:
	std::string		label_;
	PointSet*		point_set_;
	Plane3d			plane_;
	Color			color_;

	std::vector<unsigned int>	boundary_;
	std::vector<vec3>			facet_; 

	VertexGroup*			parent_;
	std::set<VertexGroup*>	children_;

	bool			visible_;
	bool			highlighted_;
};


class VertexGroupCmpDecreasing
{
public:
	VertexGroupCmpDecreasing() {}

	bool operator()(const VertexGroup* g0, const VertexGroup* g1) const {
		return g0->size() > g1->size();
	}
};


class VertexGroupCmpIncreasing
{
public:
	VertexGroupCmpIncreasing() {}

	bool operator()(const VertexGroup* g0, const VertexGroup* g1) const {
		return g0->size() < g1->size();
	}
};



inline std::vector<VertexGroup*> VertexGroup::children() {
	std::vector<VertexGroup*> tmp(children_.begin(), children_.end());
	return tmp;
}

inline void VertexGroup::set_children(const std::vector<VertexGroup*>& chld) {
	children_.clear();
	children_.insert(chld.begin(), chld.end());

	for (unsigned int i = 0; i < chld.size(); ++i)	{
		chld[i]->set_point_set(this->point_set());
		chld[i]->set_color(color());
	}
}

inline void VertexGroup::add_child(VertexGroup* g) {
	g->set_point_set(this->point_set());
	children_.insert(g);

	g->set_color(color());
}

inline void VertexGroup::remove_child(VertexGroup* g) {
	children_.erase(g);
}

inline void VertexGroup::remove_children() {
	children_.clear();
}

inline void VertexGroup::delete_children() {
	for (std::set<VertexGroup*>::iterator it = children_.begin(); it != children_.end(); ++it) {
		delete (*it);
	}
	children_.clear();
}


#endif
