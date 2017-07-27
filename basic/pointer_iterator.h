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


#ifndef _POINTER_ITERATOR_H_
#define _POINTER_ITERATOR_H_



namespace Iterators {

	//_________________________________________________________

	/**
	* PointerIterator makes it possible to iterate on a container
	* of pointers and automatically does the indirections.
	*/
	template <class CONTAINER, class T>
	class PointerIterator 
	{
	public:
		typedef PointerIterator<CONTAINER, T>    thisclass ;
		typedef typename CONTAINER::iterator     inner_iterator ;

		PointerIterator(const inner_iterator& rhs) : iterator_(rhs) {}

		thisclass& operator++() {
			iterator_++ ;
			return *this ;
		}

		thisclass& operator++(int) {
			iterator_++ ;
			return *this ;
		}

		bool operator==(const thisclass& rhs) {
			return rhs.iterator_ == iterator_ ;
		}

		bool operator!=(const thisclass& rhs) {
			return rhs.iterator_ != iterator_ ;
		}

		operator T*()   { return *iterator_; }
		T* operator->() { return *iterator_; }
		T& operator*()  { return **iterator_; }

		//______________free function____________________

		friend bool operator<(const thisclass& lhs, const thisclass& rhs) {
			return lhs.iterator_ < rhs.iterator_;
		}

	private:
		inner_iterator iterator_ ;
	} ;


	//_________________________________________________________

	/**
	* PointerConstIterator makes it possible to iterate on a container
	* of pointers and automatically does the indirections.
	*/
	template <class CONTAINER, class T> 
	class PointerConstIterator 
	{
	public:
		typedef PointerConstIterator<CONTAINER, T> thisclass ;
		typedef typename CONTAINER::const_iterator inner_iterator ;

		PointerConstIterator(
			const inner_iterator& rhs
			) : iterator_(rhs) {
		}

		thisclass& operator++() {
			iterator_++ ;
			return *this ;
		}

		thisclass& operator++(int) {
			iterator_++ ;
			return *this ;
		}

		bool operator==(const thisclass& rhs) {
			return rhs.iterator_ == iterator_ ;
		}

		bool operator!=(const thisclass& rhs) {
			return rhs.iterator_ != iterator_ ;
		}

		operator const T*() { return *iterator_; }
		const T* operator->() { return *iterator_; }
		const T& operator*() { return **iterator_; }

		//______________free function____________________

		friend bool operator<(const thisclass& lhs, const thisclass& rhs) {
			return lhs.iterator_ < rhs.iterator_;
		}

	private:
		inner_iterator iterator_ ;
	} ;

}


#endif

