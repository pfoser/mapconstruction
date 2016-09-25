package mapconstruction2;

/**
 * Frechet-based map construction 2.0 Copyright 2013 Mahmuda Ahmed and Carola Wenk
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  ------------------------------------------------------------------------
 *
 *  This software is based on the following article. Please cite this article when using this code
 * as part of a research publication:
 *
 *  Mahmuda Ahmed and Carola Wenk, "Constructing Street Networks from GPS Trajectories", European
 * Symposium on Algorithms (ESA): 60-71, Ljubljana, Slovenia, 2012
 *
 *  ------------------------------------------------------------------------
 *
 * Author: Mahmuda Ahmed Filename: Vertex.java
 *
 */


import java.util.ArrayList;
import java.util.List;

/**
 * An object that represents a point in 3D. 
 */

public class Vertex {

  private double x; // x coordinate(meters) after Mercator projection
  private double y; // y coordinate(meters) after Mercator projection
  private double z; // z coordinate(meters), it represents altitude
  private double lat; // latitude in WGS84
  private double lng; // longitude in WGS84
  private double alt; // altitude in meters

  /**
   * Contains the indices of adjacent vertices.
   */
  private List<Integer> adjacencyList;

  private boolean done = false;

  /**
   * Timestamp in milliseconds, this field is used when a pose is represented as a list of vertices.
   */
  private double timestamp = -1;

  // TODO(Mahmuda): Better to have static factory methods instead of constructor overloading.

  public Vertex() {
    this.adjacencyList = new ArrayList<Integer>();
    this.done = false;
  }

  public Vertex(double x, double y, double z) {
    this();
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public Vertex(double x, double y, double z, double timestamp) {
    this(x, y, z);
    this.timestamp = timestamp;
  }

  public Vertex(double lat, double lng, double alt, double x, double y, double z) {
    this(x, y, z);
    this.lat = lat;
    this.lng = lng;
    this.alt = alt;
  }

  public Vertex(double lat,
      double lng,
      double alt,
      double x,
      double y,
      double z,
      double timestamp) {
    this(lat, lng, alt, x, y, z);
    this.timestamp = timestamp;
  }

  public double getX() {
    return this.x;
  }

  public double getY() {
    return this.y;
  }

  public double getZ() {
    return this.z;
  }

  public double getLat() {
    return this.lat;
  }

  public double getLng() {
    return this.lng;
  }

  public double getAlt() {
    return this.alt;
  }

  public double norm(){
	return Math.sqrt(Math.pow(x, 2)+Math.pow(y, 2)+Math.pow(z, 2));  
  }
  public static double dotProd(Vertex vector1, Vertex vector2){
	  return vector1.getX()*vector2.getX()+vector1.getY()*vector2.getY()+vector1.getZ()*vector2.getZ();
  } 
  public int getDegree() {
    return this.adjacencyList.size();
  }

  public boolean getDone() {
    return this.done;
  }

  public double getTimestamp() {
    return this.timestamp;
  }


  public void setDone(boolean done) {
    this.done = done;
  }

  public List<Integer> getAdjacencyList() {
    return this.adjacencyList;
  }

  /**
   * Adds an element to its adjacency list.
   *
   * @param v is the value to be added in adjacency list
   */
  void addElementAdjList(int v) {
    for (int i = 0; i < this.getDegree(); i++) {
      if (this.adjacencyList.get(i).intValue() == v) {
        return;
      }
    }

    this.adjacencyList.add(new Integer(v));
  }

  /**
   * Returns the index of a vertex in the adjacency list.
   *
   * @param v the vertex we are looking for
   *
   * @return an int, the index of vertex k if found or -1 otherwise
   */
  public int getIndexAdjacent(int v) {
    return this.adjacencyList.indexOf(v);
  }

  /**
   * Returns the value in the adjacency list at index k
   *
   * @param k the index
   *
   * @return an int, the value at index k or -1 otherwise
   */

  public int getAdjacentElementAt(int k) {
    return this.adjacencyList.get(k).intValue();
  }

  /**
   * Set the adjacent vertex as value at index
   *
   * @param index the index to update
   *
   * @param value the new value at index
   */

  public void setAdjacentElementAt(int index, int value) {
    this.adjacencyList.remove(index);
    this.adjacencyList.add(index, value);
  }

  /**
   * Computes distance between two vertices.
   *
   * @param v2 the vertex with which we should compute distance from this vertex
   *
   * @return a double value which is the distance
   */

  public double dist(Vertex v2) {
    return Math.sqrt(Math.pow(this.x - v2.x, 2) + Math.pow(this.y - v2.y, 2));
  }

  /**
   * Resets a vertex's processing state.
   */

  public void reset() {
    done = false;
  }

  @Override
  public String toString() {
    return String.format("%f %f %f", this.x, this.y, this.z);
  }

  /**
   * @return a deep copy of this vertex
   */
  public Vertex deepCopy() {
    Vertex vertex =
        new Vertex(this.lat, this.lng, this.alt, this.x, this.y, this.z, this.timestamp);
    vertex.done = this.done;

    for (int i = 0; i < this.adjacencyList.size(); i++) {
      vertex.adjacencyList.add(this.adjacencyList.get(i).intValue());
    }
    return vertex;
  }

}
