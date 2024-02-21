Note:
If you don't want to create a visual-studio project just to compile this code into an executable, simply run the makefile and you will get the compiled executable. You can run the makefile from the terminal.

Summary:
this is a trajectory generation simulation assuming a three-link, planar manipulator. 

Each link-length can be individually defined, where:
[0] = link 1
[1] = link 2
[2] = link 3

the end-point is defined by a 3-element array where:
[0] x-cartesian coordinate
[1] y-cartesian coordinate
[2] orientation of link 3 with respect to the x-axis

the simulation 
(1) takes an initial end-point and final end-point, 
(2) computes a straight path between these two end-points, 
(3) outputs the joint angles and joint velocities to
(4) maintain a constant speed in translation and orientation between the two end-points


See "visualize_assigned.pdf" for an illustration of the joint angles as the "robot" goes from 
[3,3,0] to [5,5,0] with 
10 sampled locations in between the intial and final configuration.


