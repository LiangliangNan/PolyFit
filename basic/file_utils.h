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


#ifndef _BASIC_FILE_UTAILS_H_
#define _BASIC_FILE_UTAILS_H_

#include "basic_common.h"
#include <string>
#include <vector>

// "OpenSceneGraph - <osgDB/FileNameUtils>" has great implementation and documentation

namespace FileUtils {

	bool		BASIC_API is_file(const std::string& filename) ;
	bool		BASIC_API delete_file(const std::string& filename) ;

	bool		BASIC_API is_directory(const std::string& filename) ;   
	bool		BASIC_API delete_directory(const std::string& path) ;
	bool		BASIC_API create_directory(const std::string& path) ; // Warning: path should be absolute.

	void		BASIC_API get_directory_entries(const std::string& dir, std::vector<std::string>& entries);

	std::string BASIC_API get_current_working_directory() ;
	bool		BASIC_API set_current_working_directory(const std::string& path);

	/** Determines the home path for the current user. */
	std::string BASIC_API get_home_directory(void);

	bool		BASIC_API rename_file(const std::string& old_name, const std::string& new_name);

	time_t	BASIC_API get_time_stamp(const std::string& file_or_dir);
	std::string BASIC_API get_time_string(const std::string& file_or_dir);

	std::string BASIC_API convert_to_lower_case(const std::string& str);
	std::string BASIC_API convert_to_upper_case(const std::string& str);

	/** Gets the parent path from full name (Ex: /a/b/c.Ext => /a/b). */
	std::string BASIC_API dir_name(const std::string& file_name) ;
	/** Gets the extension without dot (Ex: /a/b/c.Ext => Ext). */
	std::string BASIC_API extension(const std::string& file_name) ;
	/** Gets the lowercase extension without dot (Ex: /a/b/c.Ext => ext). */
	std::string BASIC_API extension_in_lower_case(const std::string& filename);

	/** Gets file name with extension (Ex: /a/b/c.Ext => c.Ext). */
	std::string BASIC_API simple_name(const std::string& file_name) ;
	/** Gets file name without path and last extension (Ex: c:/file.ext1.ext2 => file.ext1; /a/b/c.Ext => c). */
	std::string BASIC_API base_name(const std::string& file_name) ;

	/** Gets file path without last extension (Ex: /a/b/c.Ext => /a/b/c ; file.ext1.ext2 => file.ext1). */
	std::string BASIC_API name_less_extension(const std::string& file_name);
	/** Gets file path without all extensions (Ex: /a/b/c.Ext => /a/b/c ; file.ext1.ext2 => file). */
	std::string BASIC_API name_less_all_extensions(const std::string& file_name);

	/**
	* Replaces extension of the given file with 'ext'. If the file name
	* does not have an extension, the given extension is appended.
	*/	
	std::string BASIC_API replace_extension(std::string const& file_name, std::string const& ext);

	/** Gets root part of a path ("/" or "C:"), or an empty string if none found. */
	std::string BASIC_API get_path_root(const std::string& path);
	/** Tests if path is absolute, as !get_path_root(path).empty(). */
	bool		BASIC_API is_absolute_path(const std::string& path);
	/** If 'to' is in a subdirectory of 'from' then this function returns the subpath, otherwise it just returns the file name.
	* The function does \b not automagically resolve paths as the system does, so be careful to give canonical paths.
	* However, the function interprets slashes ('/') ans backslashes ('\') as they were equal.
	*/
	std::string BASIC_API get_relative_path(const std::string& path);
	/** Removes .. and . dirs in a path */
	std::string BASIC_API get_absolute_path(const std::string& path);

	std::string BASIC_API convert_to_windows_style(const std::string& path);
	/** Converts back slashes (\) to forward slashes (/). */
	std::string BASIC_API convert_to_unix_style(const std::string& path);

	/** Get the path separator for the current platform. */
	char BASIC_API get_native_path_separator();
	/** Check if the path contains only the current platform's path separators. */
	bool BASIC_API is_native_style(const std::string& path);
	/** Convert the path to contain only the current platform's path separators. */
	std::string BASIC_API convert_to_native_style(const std::string& path);

	void		BASIC_API get_directory_entries(const std::string& dir, std::vector<std::string>& entries, bool recursive) ;
	void		BASIC_API get_files(const std::string& dir, std::vector<std::string>& files, bool recursive = false) ;
	void		BASIC_API get_subdirectories(const std::string& dir, std::vector<std::string>& subs, bool recursive = false) ;

	bool		BASIC_API copy_file(const std::string& original, const std::string& copy);
	bool		BASIC_API file_contains_string(const std::string& file_name, const std::string& x) ;

	void		BASIC_API read_file_to_string(const std::string& filename, std::string& data);
	void		BASIC_API write_string_to_file(const std::string& data, const std::string& filename);
	void		BASIC_API write_string_to_file(const char* data, int len, const std::string& filename);

	//////////////////////////////////////////////////////////////////////////

	// this is only for PolyFit
	std::string BASIC_API PolyFit_resource_directory();
}




#endif
