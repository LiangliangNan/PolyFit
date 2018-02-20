
PolyFit implements the hypothesis and selection based surface reconstruction method described in 
the following paper:
      --------------------------------------------------------------
      Liangliang Nan and Peter Wonka. 
      PolyFit: Polygonal Surface Reconstruction from Point Clouds. 
      ICCV 2017.
      --------------------------------------------------------------
Please consider citing the above paper if you use the code/program (or part of it). 

====================================================================================================

How to compile PolyFit from source code?
	- Download the source code of PolyFit from GitHub: https://github.com/LiangliangNan/PolyFit
	- Dependencies (only the versions specified here have been tested, while other versions should also work):
      -) Qt (version 5.8.0, 5.9.2). https://www.qt.io/
      -) CGAL (version 4.10). http://www.cgal.org/index.html
      -) boost (version 1.64. Not a direct dependence, but CGAL relies on it). http://www.boost.org/

    Besides, you need either Gurobi or lp_solve for solving the Integer Linear Programs formulated 
    by our method.
      -) Gurobi (version 7.5.1, 7.52, 64 bit). You may download it and obtain an academic license from 
	 http://www.gurobi.com/
      -) lp_solve (version 5.5, source code included in this distribution). 
	 http://lpsolve.sourceforge.net/5.5/

NOTE: Though lp_solve is included and can be used as an alternative solver, it can only solve very 
      small problems. It is too slow (may take you hours or days) and may fail to solve some 
      problems. For example, lp_solve takes 15 minutes for the data "Fig1.vg", while Gurobi takes 
      only 0.02 seconds. So don't doubt to use the Gurobi solver (academic licenses are free) that
      is much faster and reliable. 

Project file(s) for Visual Studio 2017 are provided.
The code should also be compiled by compilers under other platforms, but not tested yet.

Exectuable/binary files (may not be up-to-date).
      Pre-built binary files are available for 64-bit Windows (tested under Windows 10):
      https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html

====================================================================================================

How to run the program?
      Our algorithm consists of few major steps. To make life easier, this implementation provides a
      user interface with a few buttons (with numbered icons in an increasing order) corresponding 
      to these steps respectively. Please follow the hints on the screen to produce your results.

NOTE: This implementation incorporates a progress logger into the user interface. Thus, running 
      times should be (slightly) longer than what have been reported in our paper.	 
	
------------------
About the solvers:
------------------
	Two solvers (Gurobi and lp_solve) were linked against these binary files. However, the Gurobi 
	solver is more reliable and is always your first choice. The open source lp_solve solver can 
	solve very small problems only. It is too slow and may not guarantee to succeed. For example, 
	for the data "000-bld_02", Gurobi takes only 0.02 seconds, while lp_solve 15 minutes.  
	The dynamic library of Gurobi (version 7.51) is included in this distribution, but you may 
	still need to obtain a license from https://user.gurobi.com/download/licenses/free-academic. 
	The academic license for Gurobi is free.

====================================================================================================

Data.
      Test data can be downloaded from the project page of PolyFit:
      https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html

====================================================================================================

License.
      This program is free software; you can redistribute it and/or modify it under the terms of the
      GNU General Public License as published by the Free Software Foundation; either version 2 of 
      the License, or (at your option) any later version. The full text of the license can be found 
      in the accompanying License.txt file.

====================================================================================================

Should you have any questions, comments, or suggestions, please contact me at: 
liangliang.nan@gmail.com

Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/

July 18, 2017
Copyright (C) 2017 
