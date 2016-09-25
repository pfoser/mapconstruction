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
Filename: ReadFiles.java
 */
package mapmatchingbasics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

public class ReadFiles {

	public static ArrayList<String> getTokens(String data, String separator) {
		StringTokenizer strToken = new StringTokenizer(data, separator);
		ArrayList<String> tokens = new ArrayList<String>();
		while (strToken.hasMoreElements()) {
			tokens.add(strToken.nextToken());
		}
		return tokens;
	}

	public ArrayList<Vertex> readVertexFiles(String fileName, double XLow,
			double XHigh, double YLow, double YHigh) {
		ArrayList<Vertex> vList = new ArrayList<Vertex>();

		String str = "";
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));

			while ((str = in.readLine()) != null) {
				StringTokenizer strToken = new StringTokenizer(str, " ");
				Long id = new Long(strToken.nextToken());
				Double lat = new Double(strToken.nextToken());
				Double lon = new Double(strToken.nextToken());
				strToken.nextToken();
				strToken.nextToken();
				Double x = new Double(strToken.nextToken());
				Double y = new Double(strToken.nextToken());
				if (x >= XLow && x <= XHigh && y >= YLow && y <= YHigh)
					vList.add(new Vertex(id, lat, lon, x, y));
			}
			in.close();

		} catch (Exception e) {
			System.out.println(str + " I was here: " + e.toString());
		}
		return vList;
	}

	public HashMap<Long, Integer> getHashMap(ArrayList<Vertex> vList) {
		HashMap<Long, Integer> map = new HashMap<Long, Integer>();
		int i = 0;
		while (i < vList.size()) {
			map.put(vList.get(i).vertexIDOnFile, new Integer(i));
			i++;
		}

		return map;
	}

	public ArrayList<Vertex> removeDegreeZero(ArrayList<Vertex> vList) {
		int i = 0;
		ArrayList<Vertex> vList2 = new ArrayList<Vertex>();
		Vertex v;
		while (i < vList.size()) {
			v = vList.get(i);
			if (v.degree != 0) {
				v.degree = 0;
				vList2.add(v);
			}
			i++;
		}

		return vList2;
	}

	public ArrayList<Vertex> readEdgeFiles(String fileName,
			ArrayList<Vertex> vList, HashMap<Long, Integer> map) {
		try {

			BufferedWriter bwways = new BufferedWriter(new FileWriter(fileName
					+ "all"));
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String str;

			while ((str = in.readLine()) != null) {
				StringTokenizer strToken = new StringTokenizer(str, " ");

				strToken.nextToken();
				Long node1 = new Long(strToken.nextToken());
				Long node2 = new Long(strToken.nextToken());
				if (map.get(node1) != null && map.get(node2) != null) {
					int n1 = map.get(node1).intValue();
					int n2 = map.get(node2).intValue();

					if (n1 != n2) {
						vList.get(n1).addElementAdjList(n2);
						vList.get(n2).addElementAdjList(n1);
						bwways.write(vList.get(n1).x + " " + vList.get(n2).x
								+ " " + vList.get(n1).y + " " + vList.get(n2).y
								+ "\n");
					}
				}
			}
			in.close();
			bwways.close();

		} catch (Exception e) {
			System.out.println(e.toString());
			return vList;

		}
		return vList;

	}

	public static ArrayList<Vertex> loadBenchmarkMap(
			HashMap<Long, Integer> map, String vertexfile, String edgefile,
			boolean isDirected) {
		ArrayList<Vertex> vList = new ArrayList<Vertex>();

		try {

			BufferedReader in = new BufferedReader(new FileReader(vertexfile));
			String str;

			while ((str = in.readLine()) != null) {
				StringTokenizer strToken = new StringTokenizer(str, ",");
				Long node = new Long(strToken.nextToken());
				double x = Double.parseDouble(strToken.nextToken());
				double y = Double.parseDouble(strToken.nextToken());

				if (map.get(node) == null) {
					Vertex v = new Vertex(node, x, y);
					vList.add(v);
					int index = vList.indexOf(v);
					map.put(node, index);

				}
			}
			in.close();

		} catch (Exception e) {
			System.out.println(vertexfile + e.toString());
			System.exit(0);

		}

		try {

			BufferedReader in = new BufferedReader(new FileReader(edgefile));
			String str;

			while ((str = in.readLine()) != null) {
				if (str.trim().equals(""))
					continue;
				StringTokenizer strToken = new StringTokenizer(str, ",");
				strToken.nextToken();
				Long node1 = new Long(strToken.nextToken());
				Long node2 = new Long(strToken.nextToken());

				int index1 = -1, index2 = -1;
				if (map.get(node1) != null && map.get(node2) != null) {
					index1 = map.get(node1);
					index2 = map.get(node2);
					vList.get(index1).addElementAdjList(index2);
					if (!isDirected)
						vList.get(index2).addElementAdjList(index1);
				}
			}
			in.close();

		} catch (Exception e) {
			System.out.println(edgefile + e.toString());
			System.exit(0);

		}

		return vList;

	}

	/**
	 * Load map with defined bounding box
	 * */
	public static ArrayList<Vertex> loadBenchmarkMap(
			HashMap<Long, Integer> map, String vertexfile, String edgefile,
			boolean isDirected, double XLow, double XHigh, double YLow,
			double YHigh) {
		ArrayList<Vertex> vList = new ArrayList<Vertex>();

		try {

			BufferedReader in = new BufferedReader(new FileReader(vertexfile));
			String str;

			while ((str = in.readLine()) != null) {
				StringTokenizer strToken = new StringTokenizer(str, ",");
				Long node = new Long(strToken.nextToken());
				double x = Double.parseDouble(strToken.nextToken());
				double y = Double.parseDouble(strToken.nextToken());

				if (x >= XLow && x <= XHigh && y >= YLow && y <= YHigh) {

					if (map.get(node) == null) {
						Vertex v = new Vertex(node, x, y);
						vList.add(v);
						int index = vList.indexOf(v);
						map.put(node, index);

					}
				}

			}
			in.close();

		} catch (Exception e) {
			System.out.println(vertexfile + e.toString());
			System.exit(0);

		}

		try {

			BufferedReader in = new BufferedReader(new FileReader(edgefile));
			String str;

			while ((str = in.readLine()) != null) {
				if (str.trim().equals(""))
					continue;
				StringTokenizer strToken = new StringTokenizer(str, ",");
				strToken.nextToken();
				Long node1 = new Long(strToken.nextToken());
				Long node2 = new Long(strToken.nextToken());

				int index1 = -1, index2 = -1;
				if (map.get(node1) != null && map.get(node2) != null) {
					index1 = map.get(node1);
					index2 = map.get(node2);
					vList.get(index1).addElementAdjList(index2);
					if (!isDirected)
						vList.get(index2).addElementAdjList(index1);
				}
			}
			in.close();

		} catch (Exception e) {
			System.out.println(edgefile + e.toString());
			System.exit(0);

		}

		return vList;

	}

	public ArrayList<Vertex> readTeleAtlasVertices(String fileName,
			double XLow, double XHigh, double YLow, double YHigh) {

		ArrayList<Vertex> vList = new ArrayList<Vertex>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String str;

			while ((str = in.readLine()) != null) {
				StringTokenizer strToken = new StringTokenizer(str, " ");
				Long id = new Long(strToken.nextToken());
				Double x = new Double(strToken.nextToken());
				Double y = new Double(strToken.nextToken());
				if (x >= XLow && x <= XHigh && y >= YLow && y <= YHigh) {
					vList.add(new Vertex(id, x, y));
				}
				str = in.readLine();
				int lines = Integer.parseInt(str);
				int i = 0;
				while (i < lines) {
					in.readLine();
					i++;
				}

			}
			in.close();

		} catch (Exception e) {
			System.out.println("readTeleAtlasVertices " + e.toString());

			System.exit(0);
		}

		return vList;
	}

	public ArrayList<Vertex> readTeleAtlasGraph(String vertexFile,
			String edgeFile, double XLow, double XHigh, double YLow,
			double YHigh) {

		String str = new String();

		String vid1, vid2;
		Integer index1, index2;

		ArrayList<Vertex> vList = this.readTeleAtlasVertices(vertexFile, XLow,
				XHigh, YLow, YHigh);

		HashMap<Long, Integer> map = this.getHashMap(vList);

		try {

			BufferedReader in = new BufferedReader(new FileReader(edgeFile));

			while ((str = in.readLine()) != null) {
				try {

					StringTokenizer strToken = new StringTokenizer(str, " ");
					strToken.nextToken();// eid
					strToken.nextToken();// eid
					strToken.nextToken();// number

					vid1 = strToken.nextToken();
					vid2 = strToken.nextToken();
					index1 = map.get(new Long(vid1));
					index2 = map.get(new Long(vid2));

					if (index1 != null && index2 != null) {
						vList.get(index1.intValue()).addElementAdjList(
								index2.intValue());
						vList.get(index2.intValue()).addElementAdjList(
								index1.intValue());

					}

				} catch (Exception e) {
					System.out.println(e.toString());
				}

			}
			in.close();
		} catch (IOException e) {

			System.out.println("EXCEPTION HERE:" + e.toString());
			System.exit(0);

		}

		return vList;

	}

	public ArrayList<Vertex> loadOpenStreetMapGraph(String vertexFile,
			String edgeFile, double XLow, double XHigh, double YLow,
			double YHigh) {

		ArrayList<Vertex> vList = this.readVertexFiles(vertexFile, XLow, XHigh,
				YLow, YHigh);
		HashMap<Long, Integer> map = this.getHashMap(vList);

		vList = this.readEdgeFiles(edgeFile, vList, map);

		vList = this.removeDegreeZero(vList);
		map = this.getHashMap(vList);

		HashMap<String, Long> mapDuplicate = new HashMap<String, Long>();
		int i = 0;

		Integer index;
		while (i < vList.size()) {

			boolean contains = mapDuplicate.containsKey(vList.get(i)
					.getKeyString());

			if (!contains) {
				mapDuplicate.put(vList.get(i).getKeyString(),
						vList.get(i).vertexIDOnFile);
			} else {
				Long id = mapDuplicate.get(vList.get(i).getKeyString());
				index = map.get(id);
				map.put(vList.get(i).vertexIDOnFile, index);
			}

			i++;
		}

		vList = this.readEdgeFiles(edgeFile, vList, map);
		return vList;
	}

	public ArrayList<Vertex> loadTeleAtlasGraph(String vertexFile,
			String edgeFile, double XLow, double XHigh, double YLow,
			double YHigh) {
		ReadFiles readFiles = new ReadFiles();
		ArrayList<Vertex> vList = readFiles.readTeleAtlasGraph(vertexFile,
				edgeFile, XLow, XHigh, YLow, YHigh);
		return vList;
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

}
