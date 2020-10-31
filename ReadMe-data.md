This file describes the point cloud data used in the following paper:
```
Liangliang Nan and Peter Wonka. 
PolyFit: Polygonal Surface Reconstruction from Point Clouds. 
ICCV 2017.
```
Please consider citing the above paper if you use any of our data. 

---

### Data for PolyFit

The example data (along with additional test data and the reconstructed 3D models) are available at:
https://3d.bk.tudelft.nl/liangliang/publications/2017/polyfit/polyfit.html

For your own data, you can use my Mapple to extract planes. Here is the link to Mapple: 
https://3d.bk.tudelft.nl/liangliang/software.html    
After you load the point cloud to Mapple, go to 'Partition' menu and click 'Extract Primitives'. To visualize the planes, 
change the renderer from 'Plain' to 'Group' in the Rendering panel (at the left side of Mapple). You can save the planes as bvg (Binary Vertex Group) format. The ASCII format vg also works but slow. 

Below you will find the description of the file format and a simple example file.
 
---

### File format:

PolyFit assumes that planar segments have been extracted properly and are stored in vg (vertex group) format I have developed for my research projects. The general vg format allows you to save a point cloud followed by its segmentation information:
- For each point, its coordinates, normal, and color. 
- For each segment (a group of points representing a primitive or an object), its label, model parameters, color, and the 
    indices of the points that belong to this segment.

PolyFit handles planes only, so the model parameters are simply the plane parameters (e.g., a, b, c, and d in the plane 
equation ax + by + cz + d = 0).

Below is a details description of the ASCII vg format. The source code of PolyFit also contains an implementation for binary format. Please refer to 'point_set_serializer_vg.cpp' of the source code for more information.

```
// the coordinates of the points 
num_points: N   // N is an integer denoting the number of points
x1  y1  z1	// 3 floating point numbers
...
xN  yN  zN

// the colors of the points 
num_colors: N   // N is an integer denoting the number of colors (can be different from the number of points)
r1 g1 b1	// 3 floating point numbers
...
rN gN bN

// the normals of the points 

num_normals: N  // N is an integer denoting the number of normals (can be different from the number of points or colors)
nx1  ny1  nz1	// 3 floating point numbers
...
nxN  nyN  nzN

// now we store the segmentation information
num_groups: M   // M is an integer denoting the number of segments/primitives/objects in this point cloud

// now the information for the 1st segment/primitive/object
group_type: 0              // must be 0 standing for a plane (1 to 5 are other types of primitives)
num_group_parameters: 4    // must be 4 (planes are represented using 4 parameters) 
group_parameters: a b c d  // 4 floating point numbers (e.g., a, b, c, and d for a plane)
group_label: label         // the label (a string) of the 1st vertex group
group_color: r g b         // 3 floating point numbers denoting the color of the 1st vertex group
group_num_points: N        // N is an integer denoting the number of points in the 1st vertex group
id1 ... idN                // N integer numbers denoting the indices of the points in the 1st vertex group
num_children: 0            // a segment/primitive/object may contain subgroups, but for PolyFit this must be 0
...
group_type: 0              // here starts the last segment in the point cloud (similar to the 1st one)
num_group_parameters: 4    
group_parameters: a b c d
group_label: label
group_color: r g b
group_num_points: N
id1 ... idN
num_children: 0
```

---

### A tiny example file (comments are NOT part of the file):
```
num_points: 210882            // there are 210882 points in this point cloud
-0.06324 0.03597 0.04208 
...
-0.06449 0.03651 0.04043 
num_colors: 210882
1 0.04 0 
...
0 0.91 0
num_normals: 210882
-0.65573 -0.50319 0.56283
...
-0.56256 -0.51472 0.64698 
num_groups: 7                 // there are 7 segments in this point cloud
group_type: 0                 // the first segment is a plane
num_group_parameters: 4       // there are 4 parameters for the plane
group_parameters: -0.00494 -0.11430 0.99343 -0.03321  //a,b,c,and d in the plane equation ax+by+cz+d=0
group_label: unknown          // the label for the plane is ¡°unknown¡± (not given)
group_color: 0.1 0.6 0.2
group_num_point: 14433        // the plane consists of 14433 points
30493 ... 8798                // here are the indices of the 14433 points
num_children: 0
...                           // here are 5 segments/primitives/objects
group_type: 0                 // here comes the last segment. It is also a plane.
num_group_parameters: 4
group_parameters: 0.67223 0.25087 0.69654 -0.02011 
group_label: unknown
group_color: 0.7 0.2 0.4
group_num_point: 9155
38812 ... 140417 
num_children: 0
```

---

### Parameters. 

The parameters for most examples are as follows: fitting = 0.46, coverage = 0.27, and complexity = 0.27 (Note that weights in a 
wide range can produce the same results). Slightly different weights (fitting = 0.3, coverage = 0.4, and complexity = 0.3) are used for the sofa example in Figure 4(i), where the background (ground plane) has a much higher density than the object (sofa), 
thus the smaller data fitting weight. In case non-default parameters are used, these parameters are provided in the file names.

---

Should you have any questions, comments, or suggestions, please contact me at: 
liangliang.nan@gmail.com

**_Liangliang Nan_**

https://3d.bk.tudelft.nl/liangliang/

July 18, 2017

Copyright (C) 2017 
