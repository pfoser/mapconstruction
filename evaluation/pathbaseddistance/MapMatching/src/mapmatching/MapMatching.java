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
Filename: MapMatching.java
 */
package mapmatching;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Stack;
import java.util.StringTokenizer;

import mapmatchingbasics.*;

public class MapMatching {
	public static String city = "";

	public ArrayList<Vertex> readFile(String fileName) {

		ArrayList<Vertex> curves = new ArrayList<Vertex>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String str;
			while ((str = in.readLine()) != null) {
				StringTokenizer strToken = new StringTokenizer(str, " ");
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

	public void writeObjects(Object o, String file) {
		try {
			FileOutputStream fos = new FileOutputStream(file);
			ObjectOutputStream oos = new ObjectOutputStream(fos);

			oos.writeObject(o);

			oos.close();
		} catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public double[] computeReachabilityInterval(Vertex v1, Vertex v2,
			ArrayList<Vertex> curves, double x, double eps) {
		Line l1 = new Line(v1, v2);
		double prevt[] = new double[2];
		double t[] = l1.pIntersection(curves.get(0), eps, false);
		double interval[] = new double[2];
		boolean debug = false;
		prevt[0] = -1;
		if (debug)
			System.out.println("x = " + x + " " + v1.startindex + " "
					+ v2.startindex);
		if (x == 0)// initial situation
		{
			if (v2.startindex == 0) {
				// intervals of v1 and v2 starts in same cell
				interval[0] = v2.left;
				interval[1] = v2.right;
				return interval;
			} else {
				int cur = 1;
				prevt = l1.pIntersection(curves.get(cur), eps, false);
				while (cur <= v2.startindex) {

					t = l1.pIntersection(curves.get(cur), eps, false);

					if (t == null || t[0] > 1 || t[1] < 0)
						return null;

					prevt[0] = Math.max(t[0], prevt[0]);
					if (t[1] < prevt[0])
						return null;
					else
						prevt[1] = t[1];
					cur++;
				}
				interval[0] = v2.left;
				interval[1] = v2.right;
				return interval;
			}
		} else if (v1.startindex == v2.startindex)// In same cell
		{
			if (v1.left > v2.right)
				return null;
			interval[0] = Math.max(v2.left, v1.left);
			interval[1] = v2.right;
			return interval;

		} else {

			int cur = v1.startindex + 1;

			prevt = l1.pIntersection(curves.get(cur), eps, false);

			while (cur <= v2.startindex) {
				t = l1.pIntersection(curves.get(cur), eps, false);

				if (t == null || t[0] > 1 || t[1] < 0)
					return null;

				prevt[0] = Math.max(t[0], prevt[0]);
				// might have to check again
				if (t[1] < prevt[0])
					return null;
				else {
					prevt[1] = t[1];
				}
				cur++;

			}
			// modified 06.27.2013
			if (v2.startindex == curves.size() - 1) {
				interval[0] = curves.size() - 1;
				interval[1] = curves.size() - 1;
				return interval;
			}
			interval[0] = v2.left;
			interval[1] = v2.right;
			return interval;
		}

	}

	public ArrayList<Vertex> listStartVertices(ArrayList<Vertex> vList,
			ArrayList<Vertex> curves, double eps) {

		Vertex point, vertex2, vertex;
		ArrayList<Vertex> pq = new ArrayList<Vertex>();
		point = curves.get(0);
		int i = 0;
		boolean debug = false;
		while (i < vList.size()) {
			if (debug)
				System.out
						.println(i + " I was here before loop. "
								+ vList.get(i).toString() + " "
								+ vList.get(i).startindex + " "
								+ vList.get(i).endindex);

			this.computeNextInterval(vList.get(i), curves, 0, eps);

			if (debug)
				System.out
						.println(i + " I was here after loop. "
								+ vList.get(i).startindex + " "
								+ vList.get(i).endindex);
			i++;
		}
		i = 0;
		while (i < vList.size()) {
			vertex = vList.get(i);
			if (vertex.done)// && !vertex.done
			{

				i++;
				continue;
			}

			for (int k = 0; k < vertex.degree; k++) {
				vertex2 = vList.get(vertex.adjacencyList[k]);
				Line l1 = new Line(vertex, vertex2);
				double t[] = l1.pIntersection(point, eps, false);

				if (t != null && !(t[0] > 1 || t[1] < 0)) {

					if (!pq.contains(vertex)) {
						if (this.computeReachabilityInterval(vertex2, vertex,
								curves, 0, eps) != null)
							pq.add(vertex);

					}
					if (!pq.contains(vertex2)) {
						if (this.computeReachabilityInterval(vertex, vertex2,
								curves, 0, eps) != null)
							pq.add(vertex2);
					}
				}

			}
			i++;

		}

		return pq;
	}

	public Vertex computeNextInterval(Vertex v1, ArrayList<Vertex> curves,
			int index, double eps) {
		/**************************** startintex-----interval--------endindex *************************************/
		Line l1;
		v1.reachable = false;
		boolean debug = false;
		for (int i = index; i < curves.size() - 1; i++) {

			l1 = new Line(curves.get(i), curves.get(i + 1));

			double t[] = l1.pIntersection(v1, eps, false);

			if (t != null && !(t[0] > 1 || t[1] < 0)) {

				if (t[0] < 0)
					v1.left = i;
				else
					v1.left = i + t[0];
				v1.startindex = i;

				if (debug)
					System.out.println(i + " I was here next before. "
							+ v1.startindex + " " + t[0] + " " + t[1]);

				if (t[1] >= 1) {

					int j = i;
					while (t[1] >= 1 && j < curves.size() - 1) {
						l1 = new Line(curves.get(j), curves.get(j + 1));
						t = l1.pIntersection(v1, eps, false);
						if (debug)
							System.out.println(j + " " + t[0] + " " + t[1]);
						j++;
					}
					// made changes here
					if (j < curves.size() - 1)// t[1] < 1 &&
					{
						v1.right = j - 1 + Math.min(t[1], 1);
						v1.endindex = j;

					} else if (t[1] < 1) {
						v1.right = j + Math.min(t[1], 1);
						v1.endindex = j + 1;
					} else {
						v1.right = j + 1;
						v1.endindex = j + 1;
					}
				} else {

					v1.right = i + t[1];
					v1.endindex = i + 1;

				}
				if (debug)
					System.out.println("I was here next after. " + " "
							+ v1.left + " " + v1.right);
				return v1;
			}

		}

		v1.endindex = curves.size();
		v1.done = true;

		return v1;
	}

	public Object readObject(String fileName) {
		try {
			FileInputStream fis = new FileInputStream("t.tmp");
			ObjectInputStream ois = new ObjectInputStream(fis);

			Object obj = ois.readObject();

			ois.close();

			return obj;
		} catch (Exception e) {
			System.out.println(e.toString());
		}
		return null;
	}

	public double nextX(PriorityQueue<Vertex> pq, double x) {
		Stack<Vertex> stack = new Stack<Vertex>();
		double y = x;

		while (!pq.isEmpty() && pq.peek().left <= x) {
			stack.push(pq.poll());
		}

		if (!pq.isEmpty())
			y = pq.peek().left;
		else {
			y = Double.MAX_VALUE;
		}

		while (!stack.isEmpty()) {
			pq.add(stack.pop());
		}
		return y;
	}

	public void print(Vertex v) {
		System.out.println(v.toString() + " " + v.left + " " + v.right + " "
				+ v.startindex + " " + v.endindex);
	}

	public ArrayList<String> getFileNames(String fileName) {
		ArrayList<String> fileNames = new ArrayList<String>();

		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String str;
			while ((str = in.readLine()) != null) {
				fileNames.add(new String(str));

			}
			in.close();

		} catch (IOException e) {
			System.out.println(e.toString());
			System.exit(0);
		}

		return fileNames;

	}

	public boolean mapMatching(ArrayList<Vertex> vList,
			ArrayList<Vertex> curves, String cityname, double eps) {

		boolean debug = false, debugReturn = false;

		if (debug)
			System.out
					.println("*******************************************************************************");
		/****************************** Priority Queue *****************************************/
		Comparator<Vertex> comparator = new IntervalComparator();
		PriorityQueue<Vertex> pq = new PriorityQueue<Vertex>(1000, comparator);

		ArrayList<Vertex> startVertices = this.listStartVertices(vList, curves,
				eps);
		int i = 0;
		while (i < startVertices.size()) {
			if (debug && !startVertices.get(i).done) {
				System.out.println(startVertices.get(i).toString()
						+ startVertices.get(i).startindex + " "
						+ startVertices.get(i).endindex);

			}
			if (!startVertices.get(i).done) {
				if (startVertices.get(i).left == 0)
					startVertices.get(i).reachable = true;
				pq.add(startVertices.get(i));
			}
			i++;
		}

		if (pq.size() == 0)
			return false;
		Vertex vertex = pq.poll();
		double x = vertex.left;
		double interval[] = new double[2];
		while (vertex.right < curves.size()) {
			if (debug)
				System.out.println("\nCurrent Interval " + vertex.reachable
						+ " " + vertex.vertexID + " " + vertex.degree + " ("
						+ vertex.left + "," + vertex.right + ")");
			if (vertex.reachable && vertex.right >= curves.size() - 1)
				return true;
			for (int k = 0; k < vertex.degree; k++) {
				Vertex adjacentVertex = vList.get(vertex.adjacencyList[k]);

				if (debug && vertex.vertexID == 220)
					System.out.println("next Interval 0 "
							+ adjacentVertex.vertexID + " " + adjacentVertex.x
							+ " " + adjacentVertex.y + " "
							+ adjacentVertex.startindex + " "
							+ adjacentVertex.endindex + "("
							+ adjacentVertex.left + "," + adjacentVertex.right
							+ ")");

				if (!vertex.equals(adjacentVertex)
						&& adjacentVertex.left < adjacentVertex.right
						&& !adjacentVertex.done) {//
					if (debug)
						System.out.println("next Interval 1 "
								+ adjacentVertex.vertexID + " "
								+ adjacentVertex.startindex + " "
								+ adjacentVertex.endindex + "("
								+ adjacentVertex.left + ","
								+ adjacentVertex.right + ")");
					interval = this.computeReachabilityInterval(vertex,
							adjacentVertex, curves, x, eps);
					if (interval != null) {
						adjacentVertex.left = Math.max(interval[0],
								adjacentVertex.left);

						adjacentVertex.startindex = (int) Math
								.floor(adjacentVertex.left);
						if (interval[0] > adjacentVertex.right) {
							adjacentVertex.right = interval[1];
							adjacentVertex.endindex = (int) Math
									.ceil(adjacentVertex.right);
						}
						pq.remove(adjacentVertex);
						pq.add(adjacentVertex);
						adjacentVertex.reachable = true;
						if (debug)
							System.out.println("next Interval 2 "
									+ adjacentVertex.left + ","
									+ adjacentVertex.right);
					} else {

						if (debug)
							System.out.println("null interval");
					}
				} else if (!vertex.equals(adjacentVertex)
						&& adjacentVertex.done)// new
				{
					adjacentVertex.startindex = curves.size() - 1;
					interval = this.computeReachabilityInterval(vertex,
							adjacentVertex, curves, x, eps);

					if (interval != null && interval[0] == curves.size() - 1) {

						adjacentVertex.endindex = curves.size() - 1;
						adjacentVertex.startindex = (int) Math
								.floor(interval[0]);
						adjacentVertex.left = interval[0];
						adjacentVertex.right = interval[1];
						pq.remove(adjacentVertex);
						pq.add(adjacentVertex);
						adjacentVertex.reachable = true;
					} else {
						// do nothing
					}

				} else if (vertex.endindex == curves.size() - 1) {
					Line l1 = new Line(vertex, adjacentVertex);
					double t[] = l1.pIntersection(
							curves.get(curves.size() - 1), eps, false);
					if (t == null || t[0] > 1 || t[1] < 0) {

					} else {
						if (debugReturn)
							System.out.println("returned from case 4");
						return true;
					}
				}
			}

			pq.remove(vertex);
			if (!pq.isEmpty() && this.nextX(pq, vertex.left) < vertex.right) {
				if (debug)
					System.out.println("newly added " + !pq.isEmpty() + " "
							+ this.nextX(pq, vertex.left) + " " + vertex.right);

				vertex.left = this.nextX(pq, vertex.left);
				vertex.startindex = (int) Math.floor(vertex.left);

			}

			else {
				// x = vertex.right;
				if (debug)
					System.out.println("newly added " + !pq.isEmpty() + " "
							+ this.nextX(pq, vertex.left) + " " + vertex.right);
				vertex = this.computeNextInterval(vertex, curves,
						vertex.endindex, eps);
			}

			if (!vertex.done) {
				if (debug)
					System.out.println("Inserted again: " + vertex.vertexID
							+ " " + vertex.degree + " (" + vertex.left + ","
							+ vertex.right + ")");
				pq.add(vertex);
			}

			if (!pq.isEmpty()) {
				vertex = pq.poll();
				x = vertex.left;
				while (!vertex.reachable) {
					if (debug)
						System.out.println("Case 3: " + vertex.vertexID + " "
								+ vertex.degree + " " + vertex.reachable + " ("
								+ vertex.left + "," + vertex.right + ")");

					if (!pq.isEmpty()
							&& this.nextX(pq, vertex.left) < vertex.right) {

						if (debug)
							System.out.println("newly added in loop "
									+ !pq.isEmpty() + " "
									+ this.nextX(pq, vertex.left) + " "
									+ vertex.right);

						vertex.left = this.nextX(pq, vertex.left);
						vertex.startindex = (int) Math.floor(vertex.left);

					} else {
						if (debug)
							System.out.println("newly added in loop "
									+ !pq.isEmpty() + " "
									+ this.nextX(pq, vertex.left) + " "
									+ vertex.right);

						vertex = this.computeNextInterval(vertex, curves,
								vertex.endindex, eps);
					}
					if (!vertex.done) {
						if (debug)
							System.out.println("Inserted again in loop: "
									+ vertex.vertexID + " " + vertex.degree
									+ " (" + vertex.left + "," + vertex.right
									+ ")");
						pq.add(vertex);
					}

					if (!pq.isEmpty()) {
						vertex = pq.poll();
						x = vertex.left;
					} else {
						if (debugReturn)
							System.out.println("returned from case 3");
						return false;
					}

				}
			} else {
				if (debugReturn)
					System.out.println("returned from case 1");
				return false;
			}

			if (vertex.left > vertex.right) {

				if (debugReturn) {
					System.out.println("returned from case 2");
					System.out.println(vertex.toString() + " " + vertex.left
							+ " " + vertex.right);
				}
				return false;

			}

		}

		if (vertex.right >= curves.size()) {
			if (debugReturn)
				System.out.println("I am returned from here "
						+ vertex.toString() + " " + vertex.left + " "
						+ vertex.right + " " + vertex.reachable);
			return true;
		}

		// break;//test
		// }
		// System.out.println(startVertices.toString());
		return false;

	}

	public int min = 10000000;
	Boolean found = false;
	public String curveName;

	public void getEpsilon(ArrayList<Vertex> graph, ArrayList<Vertex> curves,
			String cityname, int start, int end) {
		boolean debug = false;
		/*
		 * if (curveName.equals("13.dat")) debug = true;
		 */

		int i;
		Boolean bool1;
		if (start >= end - 2 && !found) {
			found = true;
			for (i = start; i <= end; i++) {
				int k = 0;
				while (k < graph.size()) {
					graph.get(k).reset();
					k++;
				}
				bool1 = this.mapMatching(graph, curves, cityname, i);
				if (debug)
					System.out.println(start + ", " + end + " " + i + " = "
							+ bool1);
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
			bool1 = this
					.mapMatching(graph, curves, cityname, (end + start) / 2);
			if (debug)
				System.out.println(start + ", " + end + " " + (end + start) / 2
						+ " = " + bool1);
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

	public void pathSimilarity(ArrayList<Vertex> graph, File fin,
			String strInput, String cityname, int fileno) {
		boolean debug = false;
		ArrayList<Vertex> curves = new ArrayList<Vertex>();
		int min = 1, max = 1600;
		int count = 0;
		try {

			File folder = new File(strInput + "//");

			if (!folder.exists())
				folder.mkdirs();
			BufferedWriter bwways = new BufferedWriter(new FileWriter(strInput
					+ "//outputFrechet" + fileno + ".txt"));

			for (File file : fin.listFiles()) {

				curves = this.readFile(file.getAbsolutePath());
				this.curveName = file.getName();

				if (debug)
					System.out.println(curves.size() + " :"
							+ this.mapMatching(graph, curves, cityname, 225));
				else {

					if (curves.size() < 2) {
						System.out.println("I was here with size less than 2.");
						continue;
					}

					this.found = false;
					this.getEpsilon(graph, curves, cityname, min, max);
					if (this.min == max + 1) {
						this.found = false;
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
						+ curves.size() + " Frechet distance = " + this.min
						+ " length of path = " + dist);

				count++;
			}
			bwways.close();
		} catch (Exception e) {
			System.out.println(e.toString());
		}
	}

	public ArrayList<Vertex> removeDuplicates(ArrayList<Vertex> curves) {
		for (int i = 1; i < curves.size(); i++) {
			Vertex prev = curves.get(i - 1);
			Vertex cur = curves.get(i);

			if (prev.x == cur.x && prev.y == cur.y) {
				curves.remove(i);
				i--;
			}

		}
		return curves;
	}

}
