1) Add the two directories /source and /libraries to the current working path of the MATLAB
2) The source code is in the /source directory
3) The trajectory files should be in the form of: <x y timestamp>
4) You should firstly run the intersection_nodes_extraction.m file
   This file loads trajectories and generates intersection nodes
5) You should then run the tracebundle.m file
   This file loads trajectories and intersection nodes and generates road network vertices and edges.
   It generates two files tracebundle_vertices.txt, tracebundle_edges.txt.
6) For more information, please contact: karagior@gmail.com