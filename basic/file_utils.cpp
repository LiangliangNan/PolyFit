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


#include "file_utils.h"
#include "line_stream.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <ctime>	    // for get_time_string()

#ifdef WIN32
#include <windows.h>
#include <io.h>
#include <direct.h>   // for _mkdir
#include <sys/stat.h> // for _stat64
#include <shlobj.h>   // for SHGetFolderPathA
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <stdio.h>
#include <pwd.h>
#endif

/*
	some code are modified from or inspired by:
		"OpenSceneGraph - <osgDB/FileNameUtils>"
		"OpenSceneGraph - <osgDB/FileUtils>"
*/




namespace FileUtils {

	const char UNIX_PATH_SEPARATOR = '/';
	const char WINDOWS_PATH_SEPARATOR = '\\';

	static const char * const PATH_SEPARATORS = "/\\";
	static unsigned int PATH_SEPARATORS_LEN = 2;


	//_______________________OS-dependent functions__________________________


	bool is_file(const std::string& filename) {
#ifdef _WIN32
		struct _stat statbuf;
		if (::_stat(filename.c_str(), &statbuf) < 0)  // use '_wstat()' for Multi-Byte Character Set
			return false;

		if (!(statbuf.st_mode & _S_IFREG))
			return false;
#else // _WIN32
		struct stat statbuf;
		if (::stat(filename.c_str(), &statbuf) < 0)
			return false;

		if (!S_ISREG(statbuf.st_mode))
			return false;
#endif // _WIN32

		return true;
	}


	bool is_directory(const std::string& path) {
#ifdef _WIN32
		struct _stat statbuf;
		if (::_stat(path.c_str(), &statbuf) < 0)  // use '_wstat()' for Multi-Byte Character Set
			return false;

		if (!(statbuf.st_mode & _S_IFDIR))
			return false;
#else // _WIN32
		struct stat statbuf;
		if (::stat(path.c_str(), &statbuf) < 0)
			return false;

		if (!S_ISDIR(statbuf.st_mode))
			return false;
#endif // _WIN32

		return true;
	}


	bool create_directory(const std::string& dir) {
		if (is_directory(dir)) {
			std::cout << "directory \'" << dir << "\' already exists" << std::endl;
			return true;
		}

#ifdef WIN32
		if (::mkdir(dir.c_str()) < 0 ) { // use '_wmkdir()' for Multi-Byte Character Set
#else
		if (::mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP) != 0) {
#endif
			std::cerr << "could not mkdir" << dir << std::endl ;
			return false ;
		}

		return true ;
	}

	bool delete_directory(const std::string& path) {
		return (::rmdir(path.c_str()) == 0) ;
	}

	bool delete_file(const std::string& filename) {
		return (::unlink(filename.c_str()) == 0) ;  // you can also use "remove()"
	}


	std::string get_current_working_directory() {
		char buff[1024] ;
		return std::string(::getcwd(buff, 4096)) ;
	}


	bool set_current_working_directory(const std::string& path) {
		return (::chdir(path.c_str()) == 0) ; // // use '_wchdir()' for Multi-Byte Character Set
	}

	std::string get_home_directory() {
		char home_path[2048] = { 0 };

		if (*home_path != 0)
			return home_path;

		// TODO: Use HOME environment variable?

#ifdef _WIN32
		// SHGetFolderPathA seems to expect non-wide chars
		// http://msdn.microsoft.com/en-us/library/bb762181(VS.85).aspx
		// FIXME: Max length of home path?
		if (!SUCCEEDED(::SHGetFolderPathA(0, CSIDL_APPDATA, 0, 0, home_path)))
			std::cerr << "Cannot determine home directory" << std::endl;
#else // _WIN32
		uid_t user_id = ::geteuid();
		struct passwd* user_info = ::getpwuid(user_id);
		if (user_info == NULL || user_info->pw_dir == NULL)
			std::cerr << "Cannot determine home directory" << std::endl;
		std::strncpy(home_path, user_info->pw_dir, PATH_MAX);
#endif // _WIN32

		return home_path;
	}

	bool rename_file(const std::string& old_name, const std::string& new_name) {
		if(is_file(new_name)) {
			return false ;
		}
		return (::rename(old_name.c_str(), new_name.c_str()) == 0) ;
	}


	time_t get_time_stamp(const std::string& file_or_dir) {
		struct stat buffer;
		if (!stat(file_or_dir.c_str(), &buffer))
			return (buffer.st_mtime);
		return 0;
	}

	std::string get_time_string(const std::string& file_or_dir) {
		time_t stamp = get_time_stamp(file_or_dir);
		if (stamp != 0) {
			struct tm* timeinfo = localtime(&stamp);
			std::string tstr = asctime(timeinfo);
			//return tstr;
			return tstr.substr(0, tstr.length()-1); // discard the terminating null-character
		} else
			return "Unknown. Error occurred.";
	}


	void get_directory_entries(const std::string& dir, std::vector<std::string>& contents) {
		if (!is_directory(dir)) {
			std::cerr << "directory \'" << dir << " \' does not exist" << std::endl; 
		}

#if defined(WIN32) && !defined(__CYGWIN__)

		std::string path = dir + "/*.*";
		_finddata_t data;
		intptr_t handle = _findfirst(path.c_str(), &data);
		if (handle != -1) {
			do {
				std::string name = data.name;
				if (name != "." && name != "..") // "." and ".." seems always there
					contents.push_back(name);
			}
			while (_findnext(handle, &data) != -1);

			_findclose(handle);
		}
#else
		DIR *handle = opendir(dir.c_str());
		if (handle)
		{
			dirent *rc;
			while((rc = readdir(handle))!=NULL)
			{
				contents.push_back(rc->d_name);
			}
			closedir(handle);
		}
#endif
	}


	//_______________________OS-independent functions__________________________


	// 
	std::string convert_to_lower_case(const std::string& str)
	{
		std::string lowcase_str(str);
		for(std::string::iterator itr=lowcase_str.begin();
			itr!=lowcase_str.end();
			++itr)
		{
			*itr = tolower(*itr);
		}
		return lowcase_str;
	}
	std::string convert_to_upper_case(const std::string& str)
	{
		std::string lowcase_str(str);
		for(std::string::iterator itr=lowcase_str.begin();
			itr!=lowcase_str.end();
			++itr)
		{
			*itr = toupper(*itr);
		}
		return lowcase_str;
	}

	std::string extension(const std::string& file_name) {
		std::string::size_type dot = file_name.find_last_of('.');
		std::string::size_type slash = file_name.find_last_of(PATH_SEPARATORS);
		if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) 
			return std::string("");

		return std::string(file_name.begin() + dot + 1, file_name.end());
	}

	std::string extension_in_lower_case(const std::string& file_name)
	{
		return convert_to_lower_case(extension(file_name));
	}

	std::string base_name(const std::string& file_name) {
		std::string simpleName = simple_name(file_name);
		return name_less_extension( simpleName );
	}

	std::string dir_name(const std::string& file_name) {
		std::string::size_type slash = file_name.find_last_of(PATH_SEPARATORS);
		if (slash == std::string::npos) 
			return std::string();
		else 
			return std::string(file_name, 0, slash);
	}

	std::string simple_name(const std::string& file_name) {
		std::string::size_type slash = file_name.find_last_of(PATH_SEPARATORS);
		if (slash == std::string::npos) 
			return file_name;
		else 
			return std::string(file_name.begin() + slash + 1, file_name.end());
	}


	// strip one level of extension from the filename.
	std::string name_less_extension(const std::string& file_name)
	{
		std::string::size_type dot = file_name.find_last_of('.');
		std::string::size_type slash = file_name.find_last_of(PATH_SEPARATORS);        // Finds forward slash *or* back slash
		if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) 
			return file_name;

		return std::string(file_name.begin(), file_name.begin() + dot);
	}


	// strip all extensions from the filename.
	std::string name_less_all_extensions(const std::string& file_name) {
		// Finds start serach position: from last slash, or the begining of the string if none found
		std::string::size_type startPos = file_name.find_last_of(PATH_SEPARATORS);  // Finds forward slash *or* back slash
		if (startPos == std::string::npos) 
			startPos = 0;
		std::string::size_type dot = file_name.find_first_of('.', startPos);        // Finds *FIRST* dot from start pos
		if (dot == std::string::npos) 
			return file_name;

		return std::string(file_name.begin(), file_name.begin() + dot);
	}

	std::string replace_extension (std::string const& file_name, std::string const& ext)
	{
		std::size_t slashpos = file_name.find_last_of('/');
		if (slashpos == std::string::npos)
			slashpos = 0;

		std::size_t dotpos = file_name.find_last_of('.');
		if (dotpos == std::string::npos || dotpos < slashpos)
			return file_name + "." + ext;

		return file_name.substr(0, dotpos) + "." + ext;
	}

	std::string get_path_root(const std::string& path) {
		// Test for unix root
		if (path.empty()) 
			return "";
		if (path[0] == '/') 
			return "/";
		// Now test for Windows root
		if (path.length() < 2) 
			return "";
		if (path[1] == ':') {
			// We should check that path[0] is a letter, but as ':' is invalid in paths in other cases, that's not a problem.
			return path.substr(0, 2);
		}

		return "";
	}

	bool is_absolute_path(const std::string& path) {
#ifdef _WIN32
		// test for Windows
		// We should check that path[0] is a letter, but as ':' is invalid in paths in other cases, that's not a problem.
		return path.size() >= 2 && path[1] == ':'; 
		//return path.size() >= 2 && std::isalpha(path[0]) && path[1] == ':';
#else
		// Test for unix
		return path.size() >= 1 && path[0] == '/';
#endif   	
	}
	
	
	std::string get_relative_path(const std::string& from, const std::string& to) {
		std::cerr << "not implemented yet. Returning 'to' unchanged." << std::endl;
		return simple_name(to);
	}
	

	std::string get_absolute_path(const std::string& path) 
	{
#if defined(WIN32)  && !defined(__CYGWIN__)
		const int max_path_len = 2048;
		char retbuf[max_path_len];

		if (_fullpath(retbuf, path.c_str(), max_path_len) != 0)
			return retbuf;
		else {
			std::cerr << "invalid path. Returning 'path' unchanged." << std::endl;
			return path;
		}
#else

		char resolved_path[PATH_MAX];
		char* result = realpath(path.c_str(), resolved_path);

		if (result) return std::string(resolved_path);
		else return path;

#endif
	}


	std::string convert_to_windows_style(const std::string& path) {
		std::string new_fileName(path);

		std::string::size_type slash = 0;
		while ((slash = new_fileName.find_first_of(UNIX_PATH_SEPARATOR, slash)) != std::string::npos)
		{
			new_fileName[slash] = WINDOWS_PATH_SEPARATOR;
		}
		return new_fileName;
	}

	std::string convert_to_unix_style(const std::string& path) {
		std::string new_fileName(path);

		std::string::size_type slash = 0;
		while ((slash = new_fileName.find_first_of(WINDOWS_PATH_SEPARATOR, slash)) != std::string::npos)
		{
			new_fileName[slash] = UNIX_PATH_SEPARATOR;
		}

		return new_fileName;
	}

	char get_native_path_separator(){
#if defined(WIN32) && !defined(__CYGWIN__)
		return WINDOWS_PATH_SEPARATOR;
#else
		return UNIX_PATH_SEPARATOR;
#endif
	}

	bool is_native_style(const std::string& path){
#if defined(WIN32) && !defined(__CYGWIN__)
		return path.find(UNIX_PATH_SEPARATOR) == std::string::npos; // return true if no unix style slash exist
#else
		return path.find(WINDOWS_PATH_SEPARATOR) == std::string::npos; // return true if no windows style backslash exist
#endif
	}


	std::string convert_to_native_style(const std::string& path){
#if defined(WIN32) && !defined(__CYGWIN__)
		return convert_to_windows_style(path);
#else
		return convert_to_unix_style(path);
#endif
	}


	void get_directory_entries(
		const std::string& dir, std::vector<std::string>& result, bool recursive
		) 
	{
		get_directory_entries(dir, result) ;
		if(recursive) {
			for(unsigned int i=0; i<result.size(); i++) {
				std::string path = dir + "/" + result[i];
				if(is_directory(path)) {
					std::vector<std::string> entries;
					get_directory_entries(path, entries) ;
					for (unsigned int j=0; j<entries.size(); ++j) 
						result.push_back(result[i] + "/" + entries[j]);
				}
			}
		}
	}

	void get_files(const std::string& dir, std::vector<std::string>& result, bool recursive) {
		std::vector<std::string> entries ;
		get_directory_entries(dir, entries, recursive) ;
		for(unsigned int i=0; i<entries.size(); i++) {
			std::string name = dir + "/" + entries[i];
			if(is_file(name)) {
				result.push_back(name) ;
			}
		}
	}

	void get_subdirectories(const std::string& dir, std::vector<std::string>& result, bool recursive) {
		std::vector<std::string> entries ;
		get_directory_entries(dir, entries, recursive) ;
		for(unsigned int i=0; i<entries.size(); i++) {
			std::string name = dir + "/" + entries[i];
			if(is_directory(name)) {
				result.push_back(name) ;
			}
		}
	}

	bool copy_file(const std::string& original, const std::string& copy) {
		std::ifstream in(original.c_str());
		if (!in)
			return false;
		std::ofstream out(copy.c_str());
		LineInputStream lis(in);
		while(!lis.eof()) {
			lis.get_line();
			out << lis.current_line() << std::endl ;
		}
		return true;
	}

	bool file_contains_string(const std::string& file_name, const std::string& x) {
		std::ifstream in(file_name.c_str()) ;
		std::string buff ;
		while(in) {
			getline(in, buff) ;
			if (buff.find(x) != std::string::npos)
				return true ;
		}
		return false ;
	}


	void read_file_to_string(const std::string& filename, std::string& data) {
		std::ifstream in(filename.c_str(), std::ios::binary);
		if (in.fail()) {
			std::cerr << "Could not open file \'" << filename << "\'" << std::endl; 
			return;
		}

		in.seekg(0, std::ios::end);
		std::size_t length = in.tellg();
		in.seekg(0, std::ios::beg);
		data.resize(length);
		in.read(&(data[0]), length);
		in.close();
	}

	void write_string_to_file (const std::string& data, const std::string& filename) {
		write_string_to_file(&data[0], data.size(), filename);
	}

	void write_string_to_file (const char* data, int len, const std::string& filename) {
		std::ofstream out(filename.c_str(), std::ios::binary);
		if (out.fail()) {
			std::cerr << "Could not open file \'" << filename << "\'" << std::endl; 
			return;
		}
		out.write(data, len);
		out.close();
	}

	std::string PolyFit_resource_directory() {
		std::string dir = "resource";
		if (FileUtils::is_directory(dir)) {
			return "./" + dir;
		}	
		else if (FileUtils::is_directory("../resource"))
			return "../" + dir;
		else if (FileUtils::is_directory("../../resource"))
			return "../../" + dir;

		else if (FileUtils::is_directory("../src/resource"))
			return "../src/" + dir;
		else if (FileUtils::is_directory("../../src/resource"))
			return "../../src/" + dir;
		else if (FileUtils::is_directory("../../../src/resource"))
			return "../../../src/" + dir;

		return dir;
	}

}
