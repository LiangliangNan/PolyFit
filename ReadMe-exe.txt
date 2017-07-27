
This directory contains the executable files of PolyFit, which implements the hypothesis and selection based surface 
reconstruction method described in the following paper:

      Liangliang Nan and Peter Wonka. PolyFit: Polygonal Surface Reconstruction from Point Clouds. ICCV 2017.

These executables have been tested under 64 bit Windows 8. The source code of PolyFit is also available in the project page.

Please consider citing the above paper if you use this program. 

==============================================================================================================================

How to run the program?
      Two solvers (Gurobi and lp_solve) were linked against these binary files. However, the Gurobi solver is more reliable and 
      is always your first choice. The open source lp_solve solver can solve very small problems only. It is too slow and may 
      not guarantee to succeed. For example, for the data "000-bld_02", Gurobi takes only 0.02 seconds, while lp_solve 15 minutes.  
      The dynamic library of Gurobi (version 7.51) is included in this distribution, but you may still need to obtain a license 
      from https://user.gurobi.com/download/licenses/free-academic. The academic license for Gurobi is free.

      Our algorithm consists of few major steps. To make life easier, this implementation provides a user interface with a few 
      buttons (with numbered icons in an increasing order) corresponding to these steps respectively.

NOTE: This implementation incorporates a progress logger into the user interface. Thus, running times should be (slightly)
      longer than what have been reported in our paper.	  

==============================================================================================================================

License
      This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
      as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. The 
      full text of the license can be found in the accompanying License.txt file.

==============================================================================================================================

Should you have any questions, comments, suggestions, or you want to report issues, please contact me at: 

Liangliang Nan (liangliang.nan@gmail.com)
http://web.siat.ac.cn/~liangliang/


July 18, 2017
Copyright (C) 2017 