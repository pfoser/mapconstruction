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
Filename: HausdorffDistance.java
 */
package mapmatching;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.StringTokenizer;

import mapmatchingbasics.*;

public class HausdorffDistance {

	public boolean found;
	public double min;

	public ArrayList<Vertex> readFile(String fileName) {

		ArrayList<Vertex> curves = new ArrayList<Vertex>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String str;
			while ((str = in.readLine()) != null) {
				StringTokenizer strToken = new StringTokenizer(str);
				Double d = new Double(strToken.nextToken());
				double p1 = d.doubleValue();
				d = new Double(strToken.nextToken());
				double p2 = d.doubleValue();

				curves.add(new Vertex(p1, p2));

			}
			in.close();

		} catch (IOException e) {
			System.out.println(e.toString());
			System.exit(0);
		}

		return curves;
	}

	public boolean computeInterval(Edge e, ArrayList<Vertex> curves, int cur,
			double eps) {

		Line line = new Line(curves.get(cur - 1), curves.get(cur));
		return line.pIntersection(e, eps);
	}

	public Edge computeNextInterval(Edge e, ArrayList<Vertex> curves,
			int newstart, double eps) {
		/**************************** startintex-----interval--------endindex *************************************/

		boolean debug = false;

		boolean first = true;

		int startIndex = 0;
		double cstart = 0, vstart = 0;

		if (newstart >= curves.size()) {
			e.endIndex = curves.size();
			e.done = true;
			return e;
		}

		for (int i = newstart; i < curves.size(); i++) {
			boolean result = computeInterval(e, curves, i, eps);

			if (debug && first && result)
				System.out.println("i (next interval) = " + i + "  "
						+ (i - 1 + e.cstart) + " " + (i - 1 + e.cend) + " "
						+ first + " " + result);
			if (first && result) {
				startIndex = i - 1;
				cstart = e.cstart;
				vstart = e.vstart;

				first = false;

				if (e.cend < 1) {
					e.startIndex = startIndex;
					e.cstart = startIndex + cstart;
					e.vstart = vstart;
					e.cend = i - 1 + e.cend;
					e.endIndex = i;
					return e;
				}
			} else if (!first && result) {
				if (e.cend < 1) {
					e.startIndex = startIndex;
					e.cstart = startIndex + cstart;
					e.vstart = vstart;
					e.cend = i - 1 + e.cend;
					e.endIndex = i;
					return e;
				}
			} else if (!first && !result) {
				e.startIndex = startIndex;
				e.cstart = startIndex + cstart;
				e.vstart = vstart;
				e.cend = i - 1 + e.cend;
				e.endIndex = i;
				return e;
			}

		}

		if (first) {
			e.endIndex = curves.size();
			e.done = true;
		} else {
			e.startIndex = startIndex;
			e.cstart = startIndex + cstart;
			e.vstart = vstart;

			e.cend = curves.size() - 2 + e.cend;
			e.endIndex = curves.size() - 2;
		}

		return e;
	}

	public boolean mapMatchingHausdorff(ArrayList<Edge> graph,
			ArrayList<Vertex> curves, String cityname, double eps) {
		boolean debug = false;

		Comparator<Edge> comparator = new IntervalComparatorEdge();
		PriorityQueue<Edge> pq = new PriorityQueue<Edge>(21282, comparator);

		for (int i = 0; i < graph.size(); i++)

		{

			this.computeNextInterval(graph.get(i), curves, 1, eps);

			if (!graph.get(i).done) {

				pq.add(graph.get(i));
			}
		}
		if (pq.isEmpty())
			return false;

		Edge e = pq.poll();
		if (debug)
			System.out.println(curves.size() + "  " + e.cstart + " " + e.cend);
		double cend = e.cend;
		if (e.cstart > 0)
			return false;
		while (cend < curves.size()) {
			if (debug)
				System.out.println("interval " + e.cstart + " " + e.cend + " "
						+ cend);

			if (cend < e.cend) {
				cend = e.cend;
			}

			if (e.cend == curves.size() - 1)
				return true;

			this.computeNextInterval(e, curves, e.endIndex + 1, eps);

			if (!e.done)
				pq.add(e);

			if (pq.isEmpty()) {
				if (debug)
					System.out.println(cend + " queue empty.");
				return false;
			}

			e = pq.poll();

			if (e.cstart > cend) {
				if (debug)
					System.out.println("black interval " + e.cstart + " "
							+ cend);
				return false;
			}

		}

		return true;
	}

	public void getEpsilon(ArrayList<Edge> graph, ArrayList<Vertex> curves,
			String cityname, int start, int end) {
		int i;
		Boolean bool1;
		boolean debug = false;

		if (start >= end - 2) {
			found = true;
			for (i = start; i <= end; i++) {
				int k = 0;
				while (k < graph.size()) {
					graph.get(k).reset();
					k++;
				}
				bool1 = this.mapMatchingHausdorff(graph, curves, cityname, i);
				if (debug)
					System.out.println("here " + start + ", " + end + " " + i
							+ " = " + bool1);
				if (bool1) {
					min = i;
					return;
				}
			}
			min = i;
			return;
		} else {

			i = 0;
			while (i < graph.size()) {
				graph.get(i).reset();
				i++;
			}
			bool1 = this.mapMatchingHausdorff(graph, curves, cityname,
					(end + start) / 2);
			if (debug)
				System.out
						.println("here2 " + start + ", " + end + " "
								+ (end + start) / 2 + " = " + bool1 + " found="
								+ found);

			if (!bool1) {
				getEpsilon(graph, curves, cityname,
						(int) Math.ceil((end + start) / 2), end);
			} else {

				getEpsilon(graph, curves, cityname, start,
						(int) Math.ceil((end + start) / 2));
			}
		}

		return;
	}

	public void pathSimilarity(ArrayList<Edge> graph, File fin,
			String strInput, String cityname, int fileno) {
		boolean debug = false;
		ArrayList<Vertex> curves = new ArrayList<Vertex>();
		int min = 1, max = 1600;
		File file1 = new File(strInput);
		if (!file1.exists())
			file1.mkdirs();
		try {

			BufferedWriter bwways = new BufferedWriter(new FileWriter(strInput
					+ "//outputHausdorff" + fileno + ".txt"));
			int count = 0;
			for (File file : fin.listFiles()) {

				curves = this.readFile(file.getAbsolutePath());

				if (debug)
					System.out.println(curves.size()
							+ " :"
							+ this.mapMatchingHausdorff(graph, curves,
									cityname, 8));
				else {
					if (curves.size() < 2) {
						continue;
					}

					this.found = false;
					this.getEpsilon(graph, curves, cityname, min, max);

					if (this.min == max + 1) {
						this.found = false;
						int i = 0;
						while (i < graph.size()) {
							graph.get(i).reset();
							i++;
						}
						this.getEpsilon(graph, curves, cityname, max, max + 800);
					}
				}

				double dist = Math.sqrt(Math.pow(
						curves.get(0).x - curves.get(curves.size() - 1).x, 2)
						+ Math.pow(
								curves.get(0).y
										- curves.get(curves.size() - 1).y, 2));
				bwways.write(file.getName() + " " + curves.size() + " "
						+ this.min + " " + dist + "\n");

				System.out.println(count + "  " + file.getName() + " size = "
						+ curves.size() + " Hausdorff distance = " + this.min
						+ " length of path = " + dist);
				count++;
			}
			bwways.close();
		} catch (Exception e) {
			System.out.println(e.toString());
		}

	}

	public ArrayList<Edge> getGraphEdge(ArrayList<Vertex> vList) {
		ArrayList<Edge> graph = new ArrayList<Edge>();

		for (int i = 0; i < vList.size(); i++) {
			Vertex v1 = vList.get(i);
			for (int j = 0; j < v1.degree; j++) {
				Vertex v2 = vList.get(v1.adjacencyList[j]);
				if (!(v1.x == v2.x && v1.y == v2.y)) {
					Edge e = new Edge(v1, v2);
					graph.add(e);
				}
			}
		}

		return graph;
	}

}
