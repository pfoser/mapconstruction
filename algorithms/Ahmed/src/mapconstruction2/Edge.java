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
 * Author: Mahmuda Ahmed Filename: Edge.java
 *
 */


import java.util.ArrayList;
import java.util.List;


/**
 * An object that represents an edge from vertex1 to vertex2. 
 */

class Edge implements Comparable<Edge> {

  private Vertex vertex1; // first endpoint of edge
  private Vertex vertex2; // second endpoint of edge
  private double curveStart; // contains start point of white interval on curve
  private double curveEnd; // contains end point of white interval on curve
  private double edgeStart; // contains corresponding points on that edge
  private double edgeEnd; // contains corresponding points on that edge
  private Line line; // line from vertex v1 to v2
  /**
   * done is marked as true when an edge is done being compared with all segments of a curve.
   */
  private boolean done;
  /**
   * Contains start and end index of curves for white intervals corresponding to the edge.
   */
  private int curveStartIndex;
  private int curveEndIndex;

  /**
   * Contains a sorted list (in descending order) of positions indicating where to split this edge.
   */
  private List<Double> edgeSplitPositions;

  /**
   * For each split position the corresponding new vertex is saved in this list.
   */
  private List<Integer> edgeSplitVertices;

 // private static final FormattingLogger logger = FormattingLogger.getLogger(Edge.class);

  Edge(Vertex v1, Vertex v2) {
    this.vertex1 = v1;
    this.vertex2 = v2;
    edgeSplitPositions = new ArrayList<Double>();
    edgeSplitVertices = new ArrayList<Integer>();
    this.reset();
  }

  void reset() {
    this.curveStart = Double.MAX_VALUE;
    this.curveEnd = -1.0;
    this.edgeStart = Double.MAX_VALUE;
    this.edgeEnd = -1.0;
    this.line = new Line(vertex1, vertex2);
    this.done = false;
  }

  void set(Edge edge) {
    this.curveStart = edge.curveStart;
    this.curveEnd = edge.curveEnd;
    this.edgeStart = edge.edgeStart;
    this.edgeEnd = edge.edgeEnd;
    this.line = new Line(vertex1, vertex2);
    this.done = edge.done;
  }

  public Vertex getVertex1() {
    return this.vertex1;
  }

  public Vertex getVertex2() {
    return this.vertex2;
  }

  public Line getLine() {
    return this.line;
  }

  public boolean getDone() {
    return this.done;
  }

  public void setDone(boolean done) {
    this.done = done;
  }

  public double getCurveStart() {
    return this.curveStart;
  }

  public void setCurveStart(double cstart) {
    this.curveStart = cstart;
  }

  public double getCurveEnd() {
    return this.curveEnd;
  }

  public void setCurveEnd(double cend) {
    this.curveEnd = cend;
  }

  public double getEdgeStart() {
    return this.edgeStart;
  }

  public void setEdgeStart(double vstart) {
    this.edgeStart = vstart;
  }

  public double getEdgeEnd() {
    return this.edgeEnd;
  }

  public void setEdgeEnd(double vend) {
    this.edgeEnd = vend;
  }

  public int getCurveStartIndex() {
    return this.curveStartIndex;
  }

  public void setCurveStartIndex(int startIndex) {
    if (startIndex >= 0) {
      this.curveStartIndex = startIndex;
    } else {
     // logger.log(Level.SEVERE, "Invalid assignment of Edge.startIndex");
    }
  }

  public int getCurveEndIndex() {
    return this.curveEndIndex;
  }

  public void setCurveEndIndex(int endIndex) {
    this.curveEndIndex = endIndex;
  }

  public List<Integer> getEdgeSplitVertices() {
    return this.edgeSplitVertices;
  }

  public List<Double> getEdgeSplitPositions() {
    return this.edgeSplitPositions;
  }

  /**
   * Inserts a new split position if the list doesn't have it, otherwise return.
   *
   * @param position indicate the new split position
   * @param vertex indicates the vertex which should be inserted in this edge.
   *
   */

  public void addSplit(double position, int vertex) {
    int i = 0;
    //logger.log(Level.FINEST, "Inside updateSplits");
    for (i = 0; i < this.edgeSplitPositions.size(); i++) {
      if (this.edgeSplitPositions.get(i).doubleValue() == position) {
        return;
      } else if (this.edgeSplitPositions.get(i).doubleValue() > position) {
        this.edgeSplitPositions.add(i, position);
        this.edgeSplitVertices.add(i, vertex);
        return;
      }
    }
    this.edgeSplitPositions.add(position);
    this.edgeSplitVertices.add(vertex);
  }

  @Override
  public String toString() {
    return vertex1.toString() + vertex2.toString();
  }


  /**
   * Orders edges based on first their curveStart and then curveEnd.
   */
  @Override
  public int compareTo(Edge edge) {

    if (this.getCurveStart() < edge.getCurveStart()) {
      return -1;
    } else if (this.getCurveStart() > edge.getCurveStart()) {
      return 1;
    } else if (this.getCurveEnd() < edge.getCurveEnd()) {
      return -1;
    } else if (this.getCurveEnd() > edge.getCurveEnd()) {
      return 1;
    } else {
      return 0;
    }
  }


}
