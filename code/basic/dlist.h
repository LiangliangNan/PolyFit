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

#ifndef _DLIST_H_
#define _DLIST_H_


#include "assertions.h"

#include <vector>
#include <cstring>


// Note: if OGF_DLIST_CHECK is defined,
// some operations (such as destroy()) can
// have linear complexity (rather than constant).

#ifdef OGF_PARANOID
#define OGF_DLIST_CHECK
#endif

#ifdef OGF_DLIST_CHECK
#define ogf_dlist_assert(x) ogf_assert(x)
#else
#define ogf_dlist_assert(x)
#endif

// If your compiler barks somewhere in this file, it is
// likely that it lacks the new.h header, in that case,
// uncomment the following line:
// #define DOES_NOT_HAVE_NEW_H

#ifdef DOES_NOT_HAVE_NEW_H
inline void *operator new(size_t, void* p) { return p; }
#else
#include <new>
#endif

// The following code makes use of C++ oddities,
// enabling to separate memory management from
// object initialization and termination. 
// The correct way to call a constructor for an
// object of class T located at address addr is:
//    new(addr) T()
// which is the placement syntax for operator new.
// This means that new() uses the address addr instead
// of allocating space for the object.
// The correct way to call a destructor for an 
// object of class T located at address addr is:
//    addr->~T()
// Note that the placement syntax of operator new and
// direct calls to the destructor DO NOT allocate NOR
// deallocate memory.




//_________________________________________________________




/**
* DList is a collection and an allocator
* at the same time. It is internally represented
* by a doubly connected list. The nodes of the list
* are allocated in chunks, which both improves the
* speed and reduces the storage. 
*/

template <class T> class DList {
public:
	typedef T value_type ;
	enum { chunk_size = 128 } ;

	/**
	* Internal representation of the doubly connected
	* list. Client code does not need to access this
	* class directly.
	*/
	struct Node {
		char data_[sizeof(T)] ;
		Node* prev_ ;
		Node* next_ ;
		T* data() { return reinterpret_cast<T*>(&data_) ;  }
		const T* data() const { 
			return reinterpret_cast<const T*>(&data_) ;  
		}

		/** 
		* sanity check: allocated memory is set to zero
		* (enables uninitialized pointers to be detected
		* by VM).
		*/
		void clear_data() { Memory::clear(data_, sizeof(T)) ; }

		bool is_free() const {
			return (prev_ == nil) ;
		}

		bool is_used() const {
			return ((prev_ != nil) && (next_ != nil)) ;
		}
	} ;

protected:
	bool node_is_in_list(Node* node, Node* list) {
		for(Node* i = list->next_; i != list; i = i->next_) {
			if(i == node) {
				return true ;
			}
		}
		return false ;
	}

	void append_node_to_list(Node* node, Node* list) {
		ogf_assert(node != list) ;
		node->next_ = list ;
		list->prev_->next_ = node ;
		node->prev_ = list->prev_ ;
		list->prev_ = node ;
	}

	void remove_node_from_list(Node* node) {
		node->prev_->next_ = node->next_ ;
		node->next_->prev_ = node->prev_ ;
	}

	void init() {
		end_ = &end_node_ ;
		end_->next_ = end_ ;
		end_->prev_ = end_ ;
		inactive_end_ = &inactive_end_node_ ;
		inactive_end_->next_ = inactive_end_ ;
		inactive_end_->prev_ = inactive_end_ ;
		size_ = 0 ;
		free_list_ = nil ;
	}


public:

	class iterator {
	public:
		iterator(const iterator& rhs) : ptr_(rhs.ptr_) { }
		iterator(Node* n = nil) : ptr_(n) { }
		T& operator*()  { return *(ptr_->data()) ; }
		const T& operator*() const { return *(ptr_->data()) ; }
		T* operator->() { return ptr_->data() ; }
		const T* operator->() const { return ptr_->data() ; }
		operator T*()   { return ptr_->data() ; } 
		iterator operator++() { 
			ptr_ = ptr_->next_ ; 
			return *this ;
		}
		iterator operator++(int) {
			iterator tmp = *this ;
			++*this ;
			return tmp ;
		}
		iterator operator--() { 
			ptr_ = ptr_->prev_ ; 
			return *this ;
		}
		iterator operator--(int) {
			iterator tmp = *this ;
			--*this ;
			return tmp ;
		}

		bool operator==(const T* rhs) {
			return (ptr_ == rhs) ;
		}
		bool operator!=(const T* rhs) {
			return (ptr_ != rhs) ;
		}

		bool operator==(const iterator& rhs) {
			return (ptr_ == rhs.ptr_) ;
		}
		bool operator!=(const iterator& rhs) {
			return (ptr_ != rhs.ptr_) ;
		}

	private:
		Node* ptr_ ;
	} ;


	class const_iterator {
	public:
		const_iterator(const iterator& rhs) {
			ptr_ = reinterpret_cast<const Node*>(&*rhs) ;
		}
		const_iterator(const const_iterator& rhs) : ptr_(rhs.ptr_) { }
		const_iterator(const Node* n = nil) : ptr_(n) { }
		const T& operator*() const { return *(ptr_->data()) ; }
		const T* operator->() const { return ptr_->data() ; }
		operator const T*() const { return ptr_->data() ; } 
		const_iterator operator++() { 
			ptr_ = ptr_->next_ ; 
			return *this ;
		}
		const_iterator operator++(int) {
			const_iterator tmp = *this ;
			++*this ;
			return tmp ;
		}
		const_iterator operator--() { 
			ptr_ = ptr_->prev_ ; 
			return *this ;
		}
		const_iterator operator--(int) {
			const_iterator tmp = *this ;
			--*this ;
			return tmp ;
		}
		bool operator==(const T* rhs) {
			return (ptr_ == rhs) ;
		}
		bool operator!=(const T* rhs) {
			return (ptr_ != rhs) ;
		}
		bool operator==(const const_iterator& rhs) {
			return (ptr_ == rhs.ptr_) ;
		}
		bool operator!=(const const_iterator& rhs) {
			return (ptr_ != rhs.ptr_) ;
		}

	private:
		const Node* ptr_ ;
	} ;


	DList() {
		init() ;
	}

	/**
	* Frees all the memory associated with the data structure.
	*/
	void clear() {

		// Step 1: call the destructors.
		{for(
			Node* it = end_->next_; it != end_; it = it->next_
			) {
				// Direct call of the destructor
				// Note that this does not deallocate any
				// memory, it is just the correct way for
				// calling a destructor directly (see 
				// the comments at the beginning of this 
				// file).
				it->data()->~T() ;
		}}
		{for(
			Node* it=inactive_end_->next_; it!=inactive_end_; it=it->next_
			) {
				// Same as before
				it->data()->~T() ;
		}}

		// Step 2: deallocate memory.
		{for(unsigned int i=0; i<chunks_.size(); i++) {
			// sanity check: freed memory is set to zero
			// (enables uninitialized pointers to be detected
			// by VM).
			Memory::clear(chunks_[i], sizeof(Node) * chunk_size) ;
			delete[] (chunks_[i]) ;
		}}

		// Step 3: clear the chunks vector.
		chunks_.clear() ;

		init() ;
	}

	void clear_inactive_items() {
		Node* it = inactive_end_->next_ ;
		while(it != inactive_end_) {
			Node* next = it->next_ ;
			destroy(it->data()) ;
			it = next ;
		} 
	}

	~DList() {
		clear() ;
	}

	int size() const     { return size_ ; }
	int capacity() const { 
		return chunks_.size() * chunk_size ;
	}

	iterator begin() { return iterator(end_->next_) ; }
	iterator end() { return iterator(end_) ; }

	const_iterator begin() const { return const_iterator(end_->next_) ; }
	const_iterator end() const { return const_iterator(end_) ; }

	/*
	* creates a new element and appends it at the
	* end of the list.
	*/

	T* create() {
		if(free_list_ == nil) {
			grow() ;
		}
		Node* new_node = free_list_ ;
		free_list_ = free_list_->next_ ;

		// sanity check
		ogf_assert(new_node->is_free()) ;

		// sanity check: allocated memory is set to zero
		// (enables uninitialized pointers to be detected
		// by VM).
		new_node->clear_data() ;

		append_node_to_list(new_node, end_) ;
		size_++ ;

		// Direct call of the constructor
		// (placement syntax for operator new).
		// Note that this does not allocate any
		// memory, it is just the correct way for
		// calling a constructor directly (see 
		// the comments at the beginning of this 
		// file).
		new(new_node->data()) T() ;

		return new_node->data() ;
	}

	void destroy(T* x) {

		// Direct call of the destructor
		// Note that this does not deallocate any
		// memory, it is just the correct way for
		// calling a destructor directly (see 
		// the comments at the beginning of this 
		// file).
		x->~T() ;

		Node* node = reinterpret_cast<Node*>(x) ;

		// If this assertion fails, this
		// may be caused by a double deallocation
		// of the same element, or by a deallocation
		// of a previously deactivated element.
		ogf_assert(node->is_used()) ;
		ogf_assert(node != end_ && node != inactive_end_) ;
		remove_node_from_list(node) ;

		size_-- ;
		node->next_ = free_list_ ;
		node->prev_ = nil ;
		free_list_ = node ;

		// sanity check: deallocated memory is set to zero
		// (enables access to pointers in deallocated objects
		// to be detected by VM).
		node->clear_data() ;
	}


	/**
	* Removes an element from the list, without
	* destroying it.
	*/
	void deactivate(T* x) {
		Node* node = reinterpret_cast<Node*>(x) ;
		ogf_assert(node != end_ && node != inactive_end_) ;
		ogf_dlist_assert(node_is_in_list(node, end_)) ;
		remove_node_from_list(node) ;
		append_node_to_list(node, inactive_end_) ;
		size_-- ;
	}

	/**
	* Adds an element at the end of the list.
	* The element should have been previously
	* created using create() and removed from 
	* the list, using deactivate().
	*/
	void activate(T* x) {
		Node* node = reinterpret_cast<Node*>(x) ;
		ogf_assert(node != end_ && node != inactive_end_) ;
		ogf_dlist_assert(node_is_in_list(node, inactive_end_)) ;
		remove_node_from_list(node) ;
		append_node_to_list(node, end_) ;
		size_++ ;
	}

protected:
	void grow() {
		Node* new_chunk = new Node[chunk_size] ;

		// sanity check: allocated memory is set to zero
		// (enables uninitialized pointers to be detected
		// by VM).
		Memory::clear(new_chunk, sizeof(Node) * chunk_size) ;

		Node* last_of_chunk = &(new_chunk[chunk_size - 1]) ;
		for(Node* it = new_chunk; it != last_of_chunk; it++) {
			it->next_ = it + 1 ;
			it->prev_ = nil ;
		}
		last_of_chunk->next_ = free_list_ ;
		last_of_chunk->prev_ = nil ;
		free_list_ = new_chunk ;
		chunks_.push_back(new_chunk) ;
	}


private:
	Node end_node_ ;
	Node* end_ ;
	Node inactive_end_node_ ;
	Node* inactive_end_ ;
	int size_ ;
	Node* free_list_ ;
	std::vector<Node*> chunks_ ;

private:
	// No copy constructor nor operator= (if you really need
	//  them, I'll implement them ...
	DList(const DList<T>& rhs) ;
	DList<T>& operator=(const DList<T>& rhs) ;
} ;



#endif
