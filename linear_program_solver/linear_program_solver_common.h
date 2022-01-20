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

#ifndef SOLVER_SUITE_LINEAR_PROGRAM_SOLVER_COMMON_H
#define SOLVER_SUITE_LINEAR_PROGRAM_SOLVER_COMMON_H

// Win 32 DLL export macros
#ifdef WIN32
# ifdef SOLVER_EXPORTS
#   define SOLVER_API  __declspec(dllexport)
# else
#   define SOLVER_API  __declspec(dllimport)
# endif
# else
#define	   SOLVER_API
#endif // WIN32


#endif  // SOLVER_SUITE_LINEAR_PROGRAM_SOLVER_COMMON_H
