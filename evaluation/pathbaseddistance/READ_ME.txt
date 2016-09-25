In this implememtation of path based distance computation graph 1 will be decomposed into a set of paths and distance (Frechet/Hausdorff) will be computed with graph 2. 

Parameters
===================
(0) city name

(1) vertex file for graph 1
(2) edge file for graph 1
(3) if graph 1 is directed

(4) vertex file for graph 2
(5) edge file for graph 2
(6) if graph 2 is directed

(7) bounding box [XLow,XHigh,YLow,YHigh] (without any space) or []

(8) location where the paths will be stored
(9) location where the output files will be stored

(10) link length (LinkOne/LinkTwo/LinkThree/LineSegment)

To compute Hausdorff distance the 10th parameter is optional as it will always be computed LineSegment vs Graph.

To Compile
=============
(1) cd to the folder where src/ folder is.
For my case its: cd /home/drive1/GraphComparison/GraphDistanceBenchmark/MapMatching/MapMatching

(2) Execute the following command:
javac -d bin/ src/mapmatchingbasics/*.java src/generatepaths/*.java src/mapmatching/*.java src/benchmarkexperiments/*.java 

To Run
========

Frechet Experiment
---------------------
java -cp bin/ benchmarkexperiments.BenchmarkFrechetExperiments chicago /home/drive1/Benchmark/data/chicago/algorithms/mahmuda/chicago_mahmuda_vertices.txt /home/drive1/Benchmark/data/chicago/algorithms/mahmuda/chicago_mahmuda_edges.txt false /home/drive1/Benchmark/data/chicago/ground_truth/chicago_vertices_osm.txt /home/drive1/Benchmark/data/chicago/ground_truth/chicago_edges_osm.txt false [] /home/drive1/test/Paths/ /home/drive1/test/Results/ LinkThree

Hausdorff Experiment
------------------------
java -cp bin/ benchmarkexperiments.BenchmarkHausdorffExperiments chicago /home/drive1/Benchmark/data/chicago/algorithms/mahmuda/chicago_mahmuda_vertices.txt /home/drive1/Benchmark/data/chicago/algorithms/mahmuda/chicago_mahmuda_edges.txt false /home/drive1/Benchmark/data/chicago/ground_truth/chicago_vertices_osm.txt /home/drive1/Benchmark/data/chicago/ground_truth/chicago_edges_osm.txt false [] /home/drive1/test/Paths/ /home/drive1/test/Results/ LineSegment

Format of Output File
======================
Each output file stores following information:
(1) name of the path file
(2) complexity of the path(#of points)
(3) distance to the graph 
(4) length of the path

Process output file
=====================
To keep track of how the process is progressing the paths are stored in 5 folders (0-4) and the distance computed are stored in corresponding putput file for example: ouputFrecht0.txt contains results from the paths of folder 0. To compute the path based distance one needs to read all the output files and compute the maximum of all the distances (stored in 3rd column of output file).


