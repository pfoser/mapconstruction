/*
Path-based graph distance 1.0
Copyright 2014 Mahmuda Ahmed, K. S. Hickmann and Carola Wenk

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

------------------------------------------------------------------------

This software is based on the following article. Please cite this
article when using this code as part of a research publication:

M. Ahmed, K. S. Hickmann, and C. Wenk.
Path-based distance for street map comparison.
arXiv:1309.6131, 2013.
------------------------------------------------------------------------

Author: Mahmuda Ahmed
Filename: BenchmarkHausdorffExperments.java
 */
package benchmarkexperiments;

import generatepaths.GeneratePaths;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import mapmatching.HausdorffDistance;
import mapmatchingbasics.Edge;
import mapmatchingbasics.ReadFiles;
import mapmatchingbasics.Vertex;

public class BenchmarkHausdorffExperiments {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		GeneratePaths gp = new GeneratePaths();
		HausdorffDistance hausdorffDistance = new HausdorffDistance();

		HashMap<Long, Integer> map1 = new HashMap<Long, Integer>();
		HashMap<Long, Integer> map2 = new HashMap<Long, Integer>();

		ArrayList<Vertex> graph1 = new ArrayList<Vertex>();
		ArrayList<Vertex> graph2 = new ArrayList<Vertex>();

		String cityName = args[0];

		String vertexFile1 = args[1];
		String edgeFile1 = args[2];
		boolean isDirected1 = Boolean.parseBoolean(args[3]);

		String vertexFile2 = args[4];
		String edgeFile2 = args[5];
		boolean isDirected2 = Boolean.parseBoolean(args[6]);

		if (args[7].length() == 2) {
			graph1 = ReadFiles.loadBenchmarkMap(map1, vertexFile1, edgeFile1,
					isDirected1);
			graph2 = ReadFiles.loadBenchmarkMap(map2, vertexFile2, edgeFile2,
					isDirected2);
		} else {
			ArrayList<String> bounds = ReadFiles.getTokens(
					args[7].substring(1, args[7].length() - 1), ",");
			double XLow = Double.parseDouble(bounds.get(0));
			double XHigh = Double.parseDouble(bounds.get(1));
			double YLow = Double.parseDouble(bounds.get(2));
			double YHigh = Double.parseDouble(bounds.get(3));

			graph1 = ReadFiles.loadBenchmarkMap(map1, vertexFile1, edgeFile1,
					isDirected1, XLow, XHigh, YLow, YHigh);
			graph2 = ReadFiles.loadBenchmarkMap(map2, vertexFile2, edgeFile2,
					isDirected2, XLow, XHigh, YLow, YHigh);
		}

		
		String pathFolder = args[8];
		String resultFolder = args[9];

		String linkLength = "LineSegment";
		File pathfile = new File(pathFolder);

		if (!pathfile.exists())
			pathfile.mkdirs();

		gp.generatePathsLinkLength(graph1, pathFolder, linkLength);

		File folder = new File(pathFolder + linkLength);

		System.out.println("processing linkLength " + linkLength + " paths of "
				+ "...");
		Long start_time = System.currentTimeMillis();
		
		ArrayList<Edge> eGraph = hausdorffDistance.getGraphEdge(graph2);
		
		int count = 0;
		for (int l = 0; l < 5; l++) {
			File file2 = new File(pathFolder + folder.getName() + "//" + l);

			if (file2.exists()) {
				hausdorffDistance.pathSimilarity(eGraph, file2, resultFolder + linkLength,
						cityName, l);
				count += (file2.listFiles()).length;
			} else {
				System.out.println("Folder doesn't exits..." + pathFolder
						+ linkLength);
			}

		}
		Long end_time = System.currentTimeMillis();
		try {
			BufferedWriter bwways = new BufferedWriter(new FileWriter(
					resultFolder + "log.txt", true));

			bwways.write(cityName + " " + linkLength + " " + count + " "
					+ (end_time - start_time) / 60000.0 + "\n");

			bwways.close();
		} catch (Exception e) {
			System.out.println(e.toString());
		}

	}

}
