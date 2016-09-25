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
 * Author: Mahmuda Ahmed Filename: Line.java
 *
 */



/**
 * An object that represents a line segment between vertex p1 and p2.
 */
public class Line {

  private Vertex p1; // first endpoint
  private Vertex p2; // second endpoint
  private double xdiff; // difference between x-coordinates of p1 and p2
  private double ydiff; // difference between y-coordinates of p1 and p2
  private double zdiff; // difference between z-coordinates of p1 and p2
  private double c; // y-intersect of line equation
  private double m; // slope of line
  private double theta; // smallest angle between x-axis and this Line object
 // private static final FormattingLogger logger = FormattingLogger.getLogger(Line.class);

  /**
   * Constructor for a Line object.
   *
   * @param p1 first endpoint
   *
   * @param p2 second endpoint
   */
  public Line(Vertex p1, Vertex p2) {
    this.p1 = p1;
    this.p2 = p2;

    this.xdiff = p2.getX() - p1.getX();
    this.ydiff = p2.getY() - p1.getY();
    this.zdiff = p2.getZ() - p1.getZ();

    /*
     * y = mx + c this equation was used unless this line is parallel to y axis.
     *
     * y1 = mx1 + c y2 = mx2 + c
     *
     * c = ((y1+y2)-m(x1+x2))/2 m = (y2-y1)/(x2-x1)
     */

    if (this.xdiff != 0.0) {
      this.m = this.ydiff / this.xdiff;

      //R2Vector vector1 = new R2Vector(p1.getX(), p1.getY());
      //R2Vector vector2 = new R2Vector(p2.getX(), p2.getY());
      Vertex vector = new Vertex(p2.getX()-p1.getX(),p2.getY()-p1.getY(),0);
      double angle1 = Vertex.dotProd(vector, new Vertex(1.0, 0.0,0.0))
          / vector.norm();
      double cosAngle1 = Math.acos(angle1);
      if (this.ydiff >= 0) {
        this.theta = Math.toDegrees(cosAngle1);
      } else {
        this.theta = -1 * Math.toDegrees(cosAngle1);
      }

      this.c = ((p1.getY() + p2.getY()) - this.m * (p1.getX() + p2.getX())) / 2.0;

    } else if (this.ydiff > 0.0) {
      this.theta = Math.toDegrees(Math.PI / 2.0);
    } else {
      this.theta = -Math.toDegrees(Math.PI / 2.0);
    }
  }

  public Vertex getP1() {
    return this.p1;
  }

  public Vertex getP2() {
    return this.p2;
  }

  public double getXdiff() {
    return this.xdiff;
  }

  public double getYdiff() {
    return this.ydiff;
  }

  public double getZdiff() {
    return this.zdiff;
  }

  public double getM() {
    return this.m;
  }

  public double getC() {
    return this.c;
  }

  public double getTheta() {
    return this.theta;
  }

  public void setM(double m) {
    this.m = m;
  }

  public void setC(double c) {
    this.c = c;
  }

  public void setTheta(double theta) {
    this.theta = theta;
  }

  @Override
  public String toString() {

    return String.format("[ %f %f %f; %f %f %f]",
        this.p1.getX(),
        this.p1.getY(),
        this.p1.getZ(),
        this.p2.getX(),
        this.p2.getY(),
        this.p2.getZ());

  }

  /**
   * Compute average altitude of points p1 and p2
   */
  public double avgAltitude() {
    return (this.p1.getZ() + this.p2.getZ()) / 2.0;
  }

  /**
   * Compute intersection of eps-disc around p and line
   *
   * @param p the center of the disc
   *
   * @param eps the radius of the disc
   *
   * @return a double[2] containing two intersection points or null when they don't intersect.
   */

  public double[] pIntersection(Vertex p, double eps, boolean precisionCompromise) {

    /*
     * Line 1: x1+(x2-x1)*t = x; x1+xdiff*t = x y1+(y2-y1)*t = y; y1+ydiff*t = y
     *
     * Equation of disc around vertex p: (x-a)^2+(y-b)^2=eps^2
     * (x1+xdiff*t-a)^2+(y1+ydiff*t-b)^2=eps^2 (xdiff^2+ydiff^2)t^2 + (2(x1-a)xdiff+2(y1-b)ydiff)t +
     * (x1-a)^2+(y1-b)^2-eps^2=0
     *
     * Quadratic solution for t, gives us intersection of disc points with line
     *
     * t = (-(2(x1-a)xdiff+2(y1-b)ydiff)) +- sqrt((2(x1-a)xdiff+2(y1-b)ydiff))^2 -
     * 4(xdiff^2+ydiff^2)(x1^2+y1^2-eps^2)))/(2*(xdiff^2+ydiff^2))
     *
     * a*t^2 + b*t + c = 0
     */

    double t[] = new double[2];
    double b = 2 * ((p1.getX() - p.getX()) * this.xdiff + (p1.getY() - p.getY()) * this.ydiff);
    double c = (p1.getX() - p.getX()) * (p1.getX() - p.getX()) + (p1.getY() - p.getY())
        * (p1.getY() - p.getY()) - eps * eps;
    double a = this.xdiff * this.xdiff + this.ydiff * this.ydiff;

    if (a == 0) {
      return null;
    }

    double determinant = b * b - 4 * a * c;
    final double dist = Math.sqrt(a);
    if (determinant >= 0) {
      t[1] = (-b + Math.sqrt(determinant)) / (2 * a);
      t[0] = (-b - Math.sqrt(determinant)) / (2 * a);
    } else if (precisionCompromise) {
      double newEps = eps + 0.1;
      double newC = (p1.getX() - p.getX()) * (p1.getX() - p.getX()) + (p1.getY() - p.getY())
          * (p1.getY() - p.getY()) - newEps * newEps;
      double newDeterminant = b * b - 4 * a * newC;
      if (newDeterminant >= 0) {
        t[1] = (-b + Math.sqrt(newDeterminant)) / (2 * a);
        t[0] = (-b - Math.sqrt(newDeterminant)) / (2 * a);
      } else {
        return null;
      }
    } else {
      return null;
    }

    double min = Math.min(t[0], t[1]);
    double max = Math.max(t[0], t[1]);

    t[0] = min;
    t[1] = max;

    return t;
  }

  /**
   * Returns intersection points of this line and two lines in Line[] (neighborhood). For detail see
   * {@link package-info}
   *
   * @return an array of four double values, first two values represents interval start and end
   *         points on this line and next two corresponds to intersections with lines[0] and
   *         lines[1] respectively.
   */
  public double[] getIntervals(Line[] lines) {
    double cIntervalStart;
    double cIntervalEnd;
    double neighborhoodStart;
    double neighborhoodEnd;

    if (lines[0].getXdiff() == 0.0) {
      // when the edge is a vertical line
      cIntervalStart = this.m * lines[0].getP1().getX() + this.c;
      cIntervalEnd = this.m * lines[1].getP1().getX() + this.c;

      neighborhoodStart = (cIntervalStart - lines[0].getP1().getY()) / lines[0].getYdiff();
      neighborhoodEnd = (cIntervalEnd - lines[1].getP1().getY()) / lines[1].getYdiff();
      if (this.ydiff != 0.0) {
        cIntervalStart = (cIntervalStart - this.getP1().getY()) / this.ydiff;
        cIntervalEnd = (cIntervalEnd - this.getP1().getY()) / this.ydiff;
      } else {
        cIntervalStart = (lines[0].getP1().getX() - this.p1.getX()) / this.xdiff;
        cIntervalEnd = (lines[1].getP1().getX() - this.p1.getX()) / this.xdiff;
      }

    } else if (this.xdiff == 0.0) {
      // when the line is a vertical line
      cIntervalStart = lines[0].getM() * this.getP1().getX() + lines[0].getC();
      cIntervalEnd = lines[1].getM() * this.getP1().getX() + lines[1].getC();

      if (lines[0].getYdiff() != 0.0) {
        neighborhoodStart = (cIntervalStart - lines[0].getP1().getY()) / lines[0].getYdiff();
        neighborhoodEnd = (cIntervalEnd - lines[1].getP1().getY()) / lines[1].getYdiff();
      } else {
        neighborhoodStart = (this.p1.getX() - lines[0].getP1().getX()) / lines[0].getXdiff();
        neighborhoodEnd = (this.p1.getX() - lines[1].getP1().getX()) / lines[1].getXdiff();
      }
      cIntervalStart = (cIntervalStart - this.getP1().getY()) / this.ydiff;
      cIntervalEnd = (cIntervalEnd - this.getP1().getY()) / this.ydiff;
    } else {

      cIntervalStart = -(this.c - lines[0].getC()) / (this.m - lines[0].getM());
      cIntervalEnd = -(this.c - lines[1].getC()) / (this.m - lines[1].getM());

      neighborhoodStart = (cIntervalStart - lines[0].getP1().getX()) / lines[0].getXdiff();
      neighborhoodEnd = (cIntervalEnd - lines[1].getP1().getX()) / lines[1].getXdiff();

      cIntervalStart = (cIntervalStart - this.getP1().getX()) / this.xdiff;
      cIntervalEnd = (cIntervalEnd - this.getP1().getX()) / this.xdiff;
    }
    double temp[] = new double[4];
    temp[0] = cIntervalStart;
    temp[1] = cIntervalEnd;
    temp[2] = neighborhoodStart;
    temp[3] = neighborhoodEnd;
    return temp;
  }

  /**
   * Sets curveStart, curveEnd, edgeStart and edgeEnd on edge.
   */
  public void setIntervalOnEdge(Edge e,
      double eps,
      double cIntervalStart,
      double cIntervalEnd,
      double vstart,
      double vend) {
    // compute intersection points on the edge

    double interval[] = new double[2];

    interval[0] = Math.max(0, cIntervalStart);

    if (vstart == -1) {

      double[] in1 = e.getLine().pIntersection(this.getVertex(interval[0]), eps, true);

      if (in1 == null) {
        //logger.log(Level.SEVERE, "Problem computing Line intersection: in1.");
        throw new RuntimeException();
      }

      if (in1[0] >= 0 && in1[0] <= 1 && in1[1] >= 0 && in1[1] <= 1) {
        vstart = (in1[0] + in1[1]) / 2;
      } else if (in1[0] >= 0 && in1[0] <= 1) {
        vstart = in1[0];
      } else if (in1[1] >= 0 && in1[1] <= 1) {
        vstart = in1[1];
      } else {
        vstart = 0;
      }
    }

    interval[1] = Math.min(1, cIntervalEnd);

    if (vend == -1) {

      double[] in2 = e.getLine().pIntersection(this.getVertex(interval[1]), eps, true);

      if (in2 == null) {
        //logger.log(Level.SEVERE, "Problem computing Line intersection: in2.");
        throw new RuntimeException();
      }

      if (in2[0] >= 0 && in2[0] <= 1 && in2[1] >= 0 && in2[1] <= 1) {
        vend = (in2[0] + in2[1]) / 2;
      } else if (in2[0] >= 0 && in2[0] <= 1) {
        vend = in2[0];
      } else if (in2[1] >= 0 && in2[1] <= 1) {
        vend = in2[1];
      } else {
        vend = 1;
      }
    }

    e.setCurveStart(interval[0]);
    e.setCurveEnd(interval[1]);
    e.setEdgeStart(vstart);
    e.setEdgeEnd(vend);

  }

  /**
   * Compute intersection of eps-region around edge e and this line.
   *
   * @param e is the edge around which we would consider eps-region
   *
   * @param eps the radius of the disc
   *
   * @return true when they have intersection or false when they don't intersect.
   */
  // TODO(mahmuda): Add unit test for this method.
  public boolean pIntersection(Edge e, double eps) {

    /*
     * Line 1: x1+(x2-x1)*t = x y1+(y2-y1)*t = y
     *
     * Line 2: y = mx + c
     *
     * y1+(y2-y1)*t = (x1+(x2-x1)*t)*m + c (y2-y1)*t - (x2-x1)*t*m = x1*m + c- y1
     *
     * t = (x1*m + c - y1)/((y2-y1)-(x2-x1)*m)
     */

    Line vline = e.getLine();

    Line lines[] = Line.getEpsilonNeighborhood(vline, eps);
    double cIntervalStart;
    double cIntervalEnd;
    double neighborhoodStart;
    double neighborhoodEnd;

    double vstart = -1;
    double vend = -1;

    if (Math.abs(this.theta - vline.getTheta()) == 0
        || Math.abs(this.theta - vline.getTheta()) == 180.0) {// For parallel lines

      double t[] = this.getTParallel(vline, eps);

      if (t == null) {
        return false;
      }

      cIntervalStart = t[0];
      cIntervalEnd = t[1];

      t = vline.pIntersection(this.getVertex(cIntervalStart), eps, true);
      if (t == null) {
        return false;
      }
      neighborhoodStart = (t[0] + t[1]) / 2.0;
      t = vline.pIntersection(this.getVertex(cIntervalEnd), eps, true);
      if (t == null) {
        return false;
      }
      neighborhoodEnd = (t[0] + t[1]) / 2.0;

    } else {
      double[] temp = getIntervals(lines);
      cIntervalStart = temp[0];
      cIntervalEnd = temp[1];
      neighborhoodStart = temp[2];
      neighborhoodEnd = temp[3];
    }


    if (cIntervalStart > cIntervalEnd) {
      double temp = cIntervalStart;
      cIntervalStart = cIntervalEnd;
      cIntervalEnd = temp;

      temp = neighborhoodStart;
      neighborhoodStart = neighborhoodEnd;
      neighborhoodEnd = temp;
    }


    // computing intersection with endpoint p1

    double[] interval1 = this.pIntersection(vline.getP1(), eps, false);

    // computing intersection with endpoint p2

    double[] interval2 = this.pIntersection(vline.getP2(), eps, false);

    double minInterval1 = 0;
    double minInterval2 = 0;
    double maxInterval1 = 1;
    double maxInterval2 = 1;


    // line doesn't intersect either of the eps-disc at end points

    if (interval1 == null && interval2 == null) {

      // intersection of line and eps-neighborhood of e is non-empty
      if (cIntervalStart > 1 || cIntervalEnd < 0) {
        return false;
      }

      // intersection of line and eps-neighborhood of e is empty
      if ((neighborhoodStart > 1 && neighborhoodEnd > 1)
          || (neighborhoodStart < 0 && neighborhoodEnd < 0)) {
        return false;
      }

      // intersection needs to be in the eps-region
      if (neighborhoodStart >= 0 && neighborhoodStart <= 1 && neighborhoodEnd >= 0
          && neighborhoodEnd <= 1) {
        cIntervalStart = Math.max(0, cIntervalStart);
        cIntervalEnd = Math.min(1, cIntervalEnd);
      } else {
        return false;
      }
    }
    // line doesn't intersect with eps-disc of first endpoint

    if (interval1 != null) {
      minInterval1 = Math.min(interval1[0], interval1[1]);
      maxInterval1 = Math.max(interval1[0], interval1[1]);

      if (((neighborhoodStart > 1 && neighborhoodEnd > 1)
          || (neighborhoodStart < 0 && neighborhoodEnd < 0))
          && (minInterval1 > 1 || maxInterval1 < 0)) {
        return false;
      }
      if (neighborhoodStart < 0) {
        if (minInterval1 <= 1) {
          cIntervalStart = Math.max(cIntervalStart, minInterval1);
        }
        if (cIntervalStart == minInterval1) {
          vstart = 0;
        }
      }

      if (neighborhoodEnd < 0) {
        if (maxInterval1 >= 0) {
          cIntervalEnd = Math.min(cIntervalEnd, maxInterval1);
        }
        if (cIntervalEnd == maxInterval1) {
          vend = 0;
        }
      }
    }

    // line doesn't intersect with eps-disc of second end point

    if (interval2 != null) {

      minInterval2 = Math.min(interval2[0], interval2[1]);
      maxInterval2 = Math.max(interval2[0], interval2[1]);

      if (((neighborhoodStart > 1 && neighborhoodEnd > 1)
          || (neighborhoodStart < 0 && neighborhoodEnd < 0))
          && (minInterval2 > 1 || maxInterval2 < 0)) {
        return false;
      }

      if (neighborhoodStart > 1) {
        if (minInterval2 <= 1) {
          cIntervalStart = Math.max(cIntervalStart, minInterval2);
        }
        if (cIntervalStart == minInterval2) {
          vstart = 1;
        }
      }

      if (neighborhoodEnd > 1) {
        if (maxInterval2 >= 0) {
          cIntervalEnd = Math.min(cIntervalEnd, maxInterval2);
        }
        if (cIntervalEnd == maxInterval2) {
          vend = 1;
        }
      }
    }

    if (interval1 == null && interval2 != null) {
      if (!(neighborhoodStart >= 0 && neighborhoodStart <= 1 && neighborhoodEnd >= 0
          && neighborhoodEnd <= 1)) {
        if (neighborhoodStart >= 0 && neighborhoodStart <= 1) {
          if (cIntervalStart > 1 && minInterval2 > 1) {
            return false;
          }
          if (cIntervalStart < 0 && maxInterval2 < 0) {
            return false;
          }
        }
        if (neighborhoodEnd >= 0 && neighborhoodEnd <= 1) {
          if (cIntervalEnd > 1 && minInterval2 > 1) {
            return false;
          }
          if (cIntervalEnd < 0 && maxInterval2 < 0) {
            return false;
          }
        }
      }
    }
    if (interval1 != null && interval2 == null) {
      if (!(neighborhoodStart >= 0 && neighborhoodStart <= 1 && neighborhoodEnd >= 0
          && neighborhoodEnd <= 1)) {

        if (neighborhoodStart >= 0 && neighborhoodStart <= 1) {
          if (cIntervalStart > 1 && minInterval1 > 1) {
            return false;
          }
          if (cIntervalStart < 0 && maxInterval1 < 0) {
            return false;
          }
        }
        if (neighborhoodEnd >= 0 && neighborhoodEnd <= 1) {
          if (cIntervalEnd > 1 && minInterval1 > 1) {
            return false;
          }
          if (cIntervalEnd < 0 && maxInterval1 < 0) {
            return false;
          }
        }
      }
    }

    if (cIntervalStart > cIntervalEnd) {
      double temp = cIntervalStart;
      cIntervalStart = cIntervalEnd;
      cIntervalEnd = temp;

      temp = vend;
      vend = vstart;
      vstart = temp;
    }

    if (cIntervalStart > 1 || cIntervalEnd < 0) {
      return false;
    }

    if (Math.max(maxInterval1, maxInterval2) < 0 || Math.min(minInterval1, minInterval2) > 1) {
      return false;
    }

    this.setIntervalOnEdge(e, eps, cIntervalStart, cIntervalEnd, vstart, vend);

    return true;

  }

  /**
   * Get a Vertex on this line with parameter t.
   *
   * @param t the parameter
   *
   * @return a vertex on this line with parameter t.
   */
  public Vertex getVertex(double t) {
    return new Vertex(this.p1.getX() + this.xdiff * t, this.p1.getY() + this.ydiff * t,
        this.p1.getZ() + this.zdiff * t);
  }

  /**
   * Computes two lines which along with two eps-discs around endpoints of the line defines the
   * boundary of eps-region around the vline.
   *
   * @param vline the original segment
   *
   * @return Line[2] containing two lines
   */
  public static Line[] getEpsilonNeighborhood(Line vline, double eps) {
    // compute the equations of boundaries of eps-region around the line
    Line lines[] = new Line[2];

    double dTheta;
    if (vline.getXdiff() != 0) {
      dTheta = Math.atan(vline.getM()) + Math.PI / 2.0;
    } else if (vline.ydiff > 0.0) {
      dTheta = Math.PI / 2.0;
    } else {
      dTheta = -Math.PI / 2.0;
    }
    double dx, dy;
    dx = eps * Math.cos(dTheta);
    dy = eps * Math.sin(dTheta);

    lines[0] = new Line(
        new Vertex(vline.getP1().getX() - dx, vline.getP1().getY() - dy, vline.getP1().getZ()),
        new Vertex(vline.getP2().getX() - dx, vline.getP2().getY() - dy, vline.getP2().getZ()));
    lines[1] = new Line(
        new Vertex(vline.getP1().getX() + dx, vline.getP1().getY() + dy, vline.getP1().getZ()),
        new Vertex(vline.getP2().getX() + dx, vline.getP2().getY() + dy, vline.getP2().getZ()));

    if (lines[0].getM() != lines[1].getM()) {
      lines[0].setM(lines[1].getM());
      lines[0].setTheta(lines[1].getTheta());
    }
    return lines;
  }

  /**
   * Computes the distance between this line and a point.
   *
   * @param p the vertex from which we will compute distance
   *
   * @return a double value containing distance
   */
  public double distance(Vertex p) {

    double distance = 0;

    if (this.xdiff != 0) {
      distance =
          Math.abs(-this.m * p.getX() + p.getY() + this.c) / Math.sqrt(Math.pow(this.m, 2) + 1);
    } else {
      distance = Math.abs(this.p1.getX() - p.getX());
    }

    double t[] = this.pIntersection(p, distance, false);

    if (t == null || t[0] > 1 || t[0] < 0) {
      return Math.min(Math.sqrt((p1.getX() - p.getX()) * (p1.getX() - p.getX())
          + (p1.getY() - p.getY()) * (p1.getY() - p.getY())), Math.sqrt((p2.getX() - p.getX())
          * (p2.getX() - p.getX()) + (p2.getY() - p.getY()) * (p2.getY() - p.getY())));
    } else {
      return distance;
    }
  }

  /**
   * Computes intersections between eps-region around this line and a line segment when two lines
   * are parallel.
   *
   * @param line is a line parallel to this line
   *
   * @return a double[2] containing two intersection points or null when they don't intersect
   */

  public double[] getTParallel(Line line, double eps) {
    double t[] = new double[2];
    double newm;
    double x1, y1, x2, y2;
    if (Math.abs(line.getTheta()) == Math.PI / 2) {
      newm = 0;
      x1 = this.p1.getX();
      y1 = line.getP1().getY();
      x2 = this.p2.getX();
      y2 = line.getP2().getY();
    } else if (Math.abs(line.getTheta()) == 0) {
      newm = 0;
      x1 = line.getP1().getX();
      y1 = this.p1.getY();
      x2 = line.getP2().getX();
      y2 = this.p2.getY();
    } else {
      newm = 1 / line.getM();
      double c1 = line.getP1().getY() + newm * line.getP1().getX();
      double c2 = line.getP2().getY() + newm * line.getP2().getX();

      x1 = (c1 - this.c) / (this.m + newm);
      y1 = this.m * x1 + this.c;

      x2 = (c2 - this.c) / (this.m + newm);
      y2 = this.m * x2 + this.c;
    }

    if (Math.sqrt((line.getP1().getX() - x1) * (line.getP1().getX() - x1)
        + (line.getP1().getY() - y1) * (line.getP1().getY() - y1)) > eps) {
      return null;
    }

    double intersection1;
    double intersection2;

    if (this.xdiff != 0.0) {
      intersection1 = (x1 - this.p1.getX()) / this.xdiff;
      intersection2 = (x2 - this.p1.getX()) / this.xdiff;
    } else {
      intersection1 = (y1 - this.p1.getY()) / this.ydiff;
      intersection2 = (y2 - this.p1.getY()) / this.ydiff;
    }
    t[0] = Math.min(intersection1, intersection2);
    t[1] = Math.max(intersection1, intersection2);

    if (t[1] < 0 || t[0] > 1) {
      return null;
    }

    t[0] = Math.max(t[0], 0);
    t[1] = Math.min(t[1], 1);
    return t;
  }

  /**
   * Computes t value on this line for Vertex v, t = 0 at p1 and t = 1 at p2.
   *
   * @return a double value.
   */
  public double tValueOnLine(Vertex v) {
    if (this.xdiff == 0) {
      return (v.getY() - this.p1.getY()) / this.ydiff;
    } else {
      return (v.getX() - this.p1.getX()) / this.xdiff;
    }
  }

  /**
   * Check if the vertex, v lies on this line.
   *
   * @return boolean true, if the vertex lies on this line or false otherwise.
   */
  public boolean onLine(Vertex v) {

    if (this.xdiff == 0 && v.getX() == this.p1.getX()) {
      return true;
    } else {
      return (v.getY() - this.m * v.getX() - this.c) == 0;
    }
  }
}