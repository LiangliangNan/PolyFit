/*
*  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
*  Copyright (C) 2000-2005 INRIA - Project ALICE
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy - levy@loria.fr
*
*     Project ALICE
*     LORIA, INRIA Lorraine,
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs.
*
* As an exception to the GPL, Graphite can be linked with the following (non-GPL) libraries:
*     Qt, SuperLU, WildMagic and CGAL
*/



#ifndef _GEOM_MAP_CIRCULATORS_H_
#define _GEOM_MAP_CIRCULATORS_H_


#include "map.h"

/* usage example:

FacetHalfedgeCirculator cir(f);
	for (; !cir.end(); ++cir) {
	Point3f p = cir->vertex()->point() + translation;
	cir->vertex()->set_point(p);
}

VertexOutHalfedgeConstCirculator cir(v);
for (; !cir->end(); ++cir) {
	const Vertex* nv = cir->vertex();
	// ... 
}
*/


namespace Circulators {

	class CirculatorBase
	{
	public:
		typedef Map::Vertex*	VertexPointer ;
		typedef Map::Halfedge*	HalfedgePointer ;
		typedef Map::Facet*		FacetPointer ;

	public:
		CirculatorBase(HalfedgePointer h)
			: end_(h)
			, run_(end_)
			, steps_(0) 
		{}

		virtual ~CirculatorBase() {}

		void reset() {
			run_ = end_;
			steps_ = 0;
		}

		// Return current vertex pointed to by circulator.
		VertexPointer vertex()				{ return run_->vertex(); }
		const VertexPointer vertex() const  { return run_->vertex(); }

		// Return current halfedge pointed to by circulator.
		HalfedgePointer halfedge()			   { return run_; }
		const HalfedgePointer halfedge() const { return run_; }

		operator HalfedgePointer() { return run_; }
		operator const HalfedgePointer () const { return run_; }

		HalfedgePointer opposite() const { return run_->opposite() ; }
		HalfedgePointer next() const	 { return run_->next() ; }
		HalfedgePointer prev() const	 { return run_->prev() ; }

		bool is_border() const { 
			return run_->facet() == nil ; 
		}

		bool is_border_edge() const {
			return run_->is_border() || run_->opposite()->is_border() ; 
		}

		/// Has circulator come full circle?
		bool end() const { return run_== end_ && steps_ > 0; }

		/// Return number of steps.
		int size() const { return steps_; }

	protected:
		HalfedgePointer	 end_;
		HalfedgePointer	 run_;
		unsigned int	 steps_;
	};


	//-----------------------------------------------------------------

	class ConstCirculatorBase
	{
	public:
		typedef const Map::Vertex*     ConstVertexPointer ;
		typedef const Map::Halfedge*   ConstHalfedgePointer ;
		typedef const Map::Facet*      ConstFacetPointer ;

	public:
		ConstCirculatorBase(ConstHalfedgePointer h)
			: end_(h)
			, run_(end_)
			, steps_(0) 
		{}

		virtual ~ConstCirculatorBase() {}

		void reset() {
			run_ = end_;
			steps_ = 0;
		}

		// Return current vertex pointed to by circulator.
		ConstVertexPointer vertex() const { return run_->vertex(); }

		// Return current halfedge pointed to by circulator.
		ConstHalfedgePointer halfedge()  const { return run_; }
		operator ConstHalfedgePointer () const { return run_; }

		ConstHalfedgePointer opposite() const { return run_->opposite() ; }
		ConstHalfedgePointer next() const	  { return run_->next() ; }
		ConstHalfedgePointer prev() const	  { return run_->prev() ; }

		bool is_border() const { 
			return run_->facet() == nil ; 
		}

		bool is_border_edge() const {
			return run_->is_border() || run_->opposite()->is_border() ;
		}

		/// Has circulator come full circle?
		bool end() const { return run_== end_ && steps_ > 0; }

		/// Return number of steps.
		int size() const { return steps_; }

	protected:
		ConstHalfedgePointer	end_;
		ConstHalfedgePointer	run_;
		unsigned int			steps_;
	};



	/** 
	* Circulator to move around a facet. 
	* A circulator is similar to an iterator, which maintains the state telling us where
	* we are on a facet.
	*/
	class FacetHalfedgeCirculator : public CirculatorBase
	{
	public:
		/**
		* Construct from a face f.
		* This iterator moves around f starting from f's key halfedge.
		*/
		FacetHalfedgeCirculator(FacetPointer f) : CirculatorBase(f->halfedge())
		{}

		/** 
		* Construct from a halfedge h.
		* This iterator moves around f indicated by f starting from h.
		*/
		FacetHalfedgeCirculator(HalfedgePointer h) : CirculatorBase(h) 
		{}

		// Return the adjacent face across the current halfedge.
		FacetPointer facet()			 { return run_->opposite()->facet(); }
		const FacetPointer facet() const { return run_->opposite()->facet(); }


		// assign from non-const circulator
		FacetHalfedgeCirculator& operator=(
			const FacetHalfedgeCirculator& rhs)
		{
			end_ = rhs.end_;
			run_ = rhs.run_;
			steps_ = rhs.steps_;
			return *this;
		}

		// Equal ?
		bool operator==(
			const FacetHalfedgeCirculator& rhs) const 
		{
			return ((end_ == rhs.end_) &&
				(run_ == rhs.run_));
		}

		// Not equal ?
		bool operator!=(
			const FacetHalfedgeCirculator& rhs) const 
		{
			return !operator==(rhs);
		}

		/// Increment circulator.
		FacetHalfedgeCirculator& operator++() {
			run_=run_->next();
			++steps_;
			return *this;
		}

		/// Increment circulator.
		FacetHalfedgeCirculator& operator++(int) {
			run_=run_->next();
			++steps_;
			return *this;
		}

		/// Decrement circulator
		FacetHalfedgeCirculator& operator--() {
			run_=run_->prev();
			++steps_;
			return *this;
		}

		/// Decrement circulator
		FacetHalfedgeCirculator& operator--(int) {
			run_=run_->prev();
			++steps_;
			return *this;
		}

		/// Return a pointer to circulator self.
		FacetHalfedgeCirculator* operator->() {
			return this;
		}
	};


	//----------------------------------------------------------------------
	/** 
	* Circulator to move around a facet. 
	* A circulator is similar to an iterator, which maintains the state telling us where
	* we are on a facet.
	*/
	class FacetHalfedgeConstCirculator : public ConstCirculatorBase
	{
	public:
		/**
		* Construct from a face f.
		* This iterator moves around f starting from f's key halfedge.
		*/
		FacetHalfedgeConstCirculator(ConstFacetPointer f)
			: ConstCirculatorBase(f->halfedge())
		{}

		/** 
		* Construct from a halfedge h.
		* This iterator moves around f indicated by f starting from h.
		*/
		FacetHalfedgeConstCirculator(ConstHalfedgePointer h) 
			: ConstCirculatorBase(h) 
		{}

		// Return the adjacent face across the current halfedge.
		ConstFacetPointer facet() const {
			return run_->opposite()->facet();
		}

		// assign from non-const circulator
		FacetHalfedgeConstCirculator& operator=(
			const FacetHalfedgeConstCirculator& rhs)
		{
			end_ = rhs.end_;
			run_ = rhs.run_;
			steps_ = rhs.steps_;
			return *this;
		}

		// Equal ?
		bool operator==(
			const FacetHalfedgeConstCirculator& rhs) const 
		{
			return ((end_ == rhs.end_) && (run_ == rhs.run_));
		}

		// Not equal ?
		bool operator!=
			(const FacetHalfedgeConstCirculator& rhs) const
		{
			return !operator==(rhs);
		}

		/// Increment circulator.
		FacetHalfedgeConstCirculator& operator++() {
			run_=run_->next();
			++steps_;
			return *this;
		}

		/// Increment circulator.
		FacetHalfedgeConstCirculator& operator++(int) {
			run_=run_->next();
			++steps_;
			return *this;
		}

		/// Decrement circulator
		FacetHalfedgeConstCirculator& operator--() {
			run_=run_->prev();
			++steps_;
			return *this;
		}

		/// Decrement circulator
		FacetHalfedgeConstCirculator& operator--(int) {
			run_=run_->prev();
			++steps_;
			return *this;
		}

		/// Return a pointer to circulator self.
		FacetHalfedgeConstCirculator* operator->() {
			return this;
		}
	};


	//---------------------------------------------------------------------

	/** 
	* Circulator for moving around a vertex.
	* This circulator makes it easy to visit all faces, edges, and vertices
	* adjacent to a given vertex.
	* NOTE:: results are halfedges radiating from the vertex! 
	*/
	class VertexOutHalfedgeCirculator : public CirculatorBase
	{
	public:
		/** 
		* Construct a vertex circulator from a vertex.
		* This constructor creates a vertex circulator which starts at a random
		* point in the one ring.
		*/
		VertexOutHalfedgeCirculator(VertexPointer v) 
			: CirculatorBase(v->halfedge()->opposite()) 
		{}

		/// Get the face of the current halfedge.
		FacetPointer facet() { 
			return run_->facet();
		}

		const FacetPointer facet() const {
			return run_->facet();
		}

		// assign from non-const circulator
		VertexOutHalfedgeCirculator& operator=(
			const VertexOutHalfedgeCirculator& rhs) 
		{
			end_ = rhs.end_;
			run_ = rhs.run_;
			steps_ = rhs.steps_;
			return *this;
		}

		// Equal ?
		bool operator==(
			const VertexOutHalfedgeCirculator& rhs) const
		{
			return ((end_ == rhs.end_) && (run_ == rhs.run_));
		}

		/// Not equal ?
		bool operator!=(
			const VertexOutHalfedgeCirculator& rhs) const
		{
			return !operator==(rhs);
		}

		// Increment circulator (ccw)
		VertexOutHalfedgeCirculator& operator++() {
			run_=run_->prev()->opposite();
			++steps_;
			return *this;
		}

		// Increment circulator (ccw)
		VertexOutHalfedgeCirculator& operator++(int) {
			run_=run_->prev()->opposite();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexOutHalfedgeCirculator& operator--() {
			run_=run_->opposite()->next();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexOutHalfedgeCirculator& operator--(int) {
			run_=run_->opposite()->next();
			++steps_;
			return *this;
		}	

		// Return a pointer to circulator self.
		VertexOutHalfedgeCirculator* operator->() {
			return this;
		}
	};


	//---------------------------------------------------------------------

	/** 
	* Circulator for moving around a vertex.
	* This circulator makes it easy to visit all faces, edges, and vertices
	* adjacent to a given vertex.
	* NOTE:: results are halfedges radiating from the vertex! 
	*/
	class VertexOutHalfedgeConstCirculator : public ConstCirculatorBase
	{
	public:
		/** 
		* Construct a vertex circulator from a vertex.
		* This constructor creates a vertex circulator which starts at a random
		* point in the one ring.
		*/
		VertexOutHalfedgeConstCirculator(ConstVertexPointer v) 
			: ConstCirculatorBase(v->halfedge()->opposite()) 
		{}

		/// Get the face of the current halfedge.
		ConstFacetPointer facet() const {
			return run_->facet(); 
		}

		// assign from non-const circulator
		VertexOutHalfedgeConstCirculator& operator=(
			const VertexOutHalfedgeConstCirculator& rhs) 
		{
			end_ = rhs.end_;
			run_ = rhs.run_;
			steps_ = rhs.steps_;
			return *this;
		}

		// Equal ?
		bool operator==(
			const VertexOutHalfedgeConstCirculator& rhs) const 
		{
			return ((end_ == rhs.end_) && (run_ == rhs.run_));
		}

		/// Not equal ?
		bool operator!=(
			const VertexOutHalfedgeConstCirculator& rhs) const 
		{
			return !operator==(rhs);
		}

		// Increment circulator (ccw)
		VertexOutHalfedgeConstCirculator& operator++() {
			run_=run_->prev()->opposite();
			++steps_;
			return *this;
		}

		// Increment circulator (ccw)
		VertexOutHalfedgeConstCirculator& operator++(int) {
			run_=run_->prev()->opposite();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexOutHalfedgeConstCirculator& operator--() {
			run_=run_->opposite()->next();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexOutHalfedgeConstCirculator& operator--(int) {
			run_=run_->opposite()->next();
			++steps_;
			return *this;
		}	

		// Return a pointer to circulator self.
		VertexOutHalfedgeConstCirculator* operator->() {
			return this;
		}
	};


	//---------------------------------------------------------------------


	/** 
	* Circulator for moving around a vertex.
	* This circulator makes it easy to visit all faces, edges, and vertices
	* adjacent to a given vertex.
	* NOTE:: results are halfedges pointing to the vertex! 
	*/
	class VertexInHalfedgeCirculator : public CirculatorBase
	{
	public:
		/** 
		* Construct a vertex circulator from a vertex.
		* This constructor creates a vertex circulator which starts at a random
		* point in the one ring.
		*/
		VertexInHalfedgeCirculator(VertexPointer v)
			: CirculatorBase(v->halfedge()) 
		{}

		// Get the face of the current halfedge.
		FacetPointer facet() {
			return run_->facet(); 
		}

		const FacetPointer facet() const {
			return run_->facet();
		}

		// assign from non-const circulator
		VertexInHalfedgeCirculator& operator=(
			const VertexInHalfedgeCirculator& rhs) {
				end_ = rhs.end_;
				run_ = rhs.run_;
				steps_ = rhs.steps_;
				return *this;
		}

		// Equal ?
		bool operator==(
			const VertexInHalfedgeCirculator& rhs) const
		{
			return ((end_ == rhs.end_) && (run_ == rhs.run_));
		}

		// Not equal ?
		bool operator!=(
			const VertexInHalfedgeCirculator& rhs) const
		{
			return !operator==(rhs);
		}


		// Increment circulator (ccw)
		VertexInHalfedgeCirculator& operator++() {
			run_=run_->opposite()->prev();
			++steps_;
			return *this;
		}

		// Increment circulator (ccw)
		VertexInHalfedgeCirculator& operator++(int) {
			run_=run_->opposite()->prev();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexInHalfedgeCirculator& operator--() {
			run_=run_->next()->opposite();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexInHalfedgeCirculator& operator--(int) {
			run_=run_->next()->opposite();
			++steps_;
			return *this;
		}	

		// Return a pointer to circulator self.
		VertexInHalfedgeCirculator* operator->() {
			return this;
		}


	};


	//---------------------------------------------------------------------


	/** 
	* Circulator for moving around a vertex.
	* This circulator makes it easy to visit all faces, edges, and vertices
	* adjacent to a given vertex.
	* NOTE:: results are halfedges pointing to the vertex! 
	*/
	class VertexInHalfedgeConstCirculator : public ConstCirculatorBase
	{
	public:
		/** 
		* Construct a vertex circulator from a vertex.
		* This constructor creates a vertex circulator which starts at a random
		* point in the one ring.
		*/
		VertexInHalfedgeConstCirculator(ConstVertexPointer v) 
			: ConstCirculatorBase(v->halfedge()) 
		{}


		// Get the face of the current halfedge.
		ConstFacetPointer facet() const {
			return run_->facet(); 
		}

		// assign from non-const circulator
		VertexInHalfedgeConstCirculator& operator=(
			const VertexInHalfedgeConstCirculator& rhs) 
		{
			end_ = rhs.end_;
			run_ = rhs.run_;
			steps_ = rhs.steps_;
			return *this;
		}

		// Equal ?
		bool operator==(
			const VertexInHalfedgeConstCirculator& rhs) const 
		{
			return ((end_ == rhs.end_) && (run_ == rhs.run_));
		}

		// Not equal ?
		bool operator!=(
			const VertexInHalfedgeConstCirculator& rhs) const
		{
			return !operator==(rhs);
		}


		// Increment circulator (ccw)
		VertexInHalfedgeConstCirculator& operator++() {
			run_=run_->opposite()->prev();
			++steps_;
			return *this;
		}

		// Increment circulator (ccw)
		VertexInHalfedgeConstCirculator& operator++(int) {
			run_=run_->opposite()->prev();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexInHalfedgeConstCirculator& operator--() {
			run_=run_->next()->opposite();
			++steps_;
			return *this;
		}

		// Decrement circulator (cw)
		VertexInHalfedgeConstCirculator& operator--(int) {
			run_=run_->next()->opposite();
			++steps_;
			return *this;
		}	

		// Return a pointer to the circulator self.
		VertexInHalfedgeConstCirculator* operator->() {
			return this;
		}

	};

}

typedef Circulators::FacetHalfedgeCirculator			FacetHalfedgeCirculator;
typedef Circulators::VertexInHalfedgeCirculator			VertexInHalfedgeCirculator;
typedef Circulators::VertexOutHalfedgeCirculator		VertexOutHalfedgeCirculator;

typedef Circulators::FacetHalfedgeConstCirculator		FacetHalfedgeConstCirculator;
typedef Circulators::VertexInHalfedgeConstCirculator	VertexInHalfedgeConstCirculator;
typedef Circulators::VertexOutHalfedgeConstCirculator	VertexOutHalfedgeConstCirculator;



#endif
