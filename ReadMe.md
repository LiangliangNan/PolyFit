PolyFit implements the hypothesis and selection based surface reconstruction method described in the following paper:
```
Liangliang Nan and Peter Wonka. 
PolyFit: Polygonal Surface Reconstruction from Point Clouds. 
ICCV 2017.
```
Please consider citing the above paper if you use the code/program (or part of it). 

=========================================================================

### Run PolyFit
- Download executable/binary files (tested on 64-bit Windows 10 only): 
  https://github.com/LiangliangNan/PolyFit/releases

  Note: The release available might not be the latest stable one. I recommend building PolyFit from the source code (see below).

- Follow the hints on the screen to produce your results.
  
  Our algorithm consists of few major steps. This demo version provides a user interface with a few buttons (with numbered icons) and screen hints corresponding to these steps.

### About the solvers
Two solvers, namely Gurobi and lp_solve, are available in PolyFit. The Gurobi solver is more reliable and is always your first choice. The open source lp_solve solver can only solve tiny problems. It is too slow and may not guarantee to succeed. For example the data "000-bld_02", Gurobi takes only 0.02 seconds, while lp_solve 15 minutes. For your convenience, the dynamic library of Gurobi is included in this distribution, but you may still need to obtain a license (free for academic use) from 
https://user.gurobi.com/download/licenses/free-academic. 
      
### About the timing
This implementation incorporates a progress logger into the user interface. Thus, running times should be (slightly) longer than what has been reported in our paper.     

=========================================================================

### Build PolyFit from source code
* Download the source code of PolyFit from GitHub: https://github.com/LiangliangNan/PolyFit
* Install the dependencies (though versions specified have been tested, other versions should also work):
  - Qt (v5.8.0, v5.9.2, v5.10.1): 
    https://www.qt.io/
  - CGAL (v4.10, v4.11.1):
    http://www.cgal.org/index.html
  - boost (v1.64. Not a direct dependence, but CGAL relies on it):
    http://www.boost.org/
  - Gurobi (v7.5.1, v7.52, 64 bit). You may download it and obtain an academic license from 
    http://www.gurobi.com/
* Build PolyFit. 
  Project files for the following IDEs (Integrated Development Environment) are provided:
  - QtCreator (for macOS, Linux, Windows, etc., but only tested on macOS): PolyFit.pro
  - Visual Studio 2017 (for Windows, only tested on 64bits Windows 10): PolyFit.sln
  
You should be able to build PolyFit on most platforms with little effort in editing the project files.

=========================================================================

### Data
Test data can be downloaded from the project page of PolyFit:
https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html

#### Data format
Please have a look at:
https://github.com/LiangliangNan/PolyFit/blob/master/ReadMe-data.txt

#### Plane extraction
PolyFit assumes that the planar segments are provided as input. 
Extracting planes has some randomness (due to the nature of RANSAC) and the data quality can vary a lot (it should be fine if some regions of the planes are missing), so I isolated this part from this demo version implementation. You can use my Mapple to extract planes from point clouds. Here is the link to Mapple: https://3d.bk.tudelft.nl/liangliang/software.html    

After you load the point cloud to Mapple, go to 'Partition' menu and click 'Extract Primitives'. To visualize the planes, change the renderer from 'Plain' to 'Group' in the Rendering panel (at the left side of Mapple). You can save the planes as bvg (**B**inary **V**ertex **G**roup) format. The ASCII format vg also works but slow.

=========================================================================

### License
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License or (at your option) any later version. The full text of the license can be found in the accompanying License.txt file.

=========================================================================

Should you have any questions, comments, or suggestions, please contact me at: 
liangliang.nan@gmail.com

**_Liangliang Nan_**

https://3d.bk.tudelft.nl/liangliang/

July 18, 2017

Copyright (C) 2017 
