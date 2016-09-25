1) Add the two directories /source and /libraries to the current working path of the MATLAB
2) The source code is in the /source directory
3) As an example, the /source directory also contains the ground truth and generated road networks as explained below:
   Gground truth data: vertices_original_osm.txt, edges_original_osm.txt
   Generated road network data: tracebundle_vertices.txt, tracebundle_edges.txt
4) You should run shortest_paths_comparison.m file
   This file generates uniformly distributed shortest paths and computes the geometric similarity/difference
   of the ground truth and generated road networks
5) For more information, please contact: karagior@gmail.com