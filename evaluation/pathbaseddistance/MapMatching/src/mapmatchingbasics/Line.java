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
Filename: Line.java
 */
package mapmatchingbasics;

public class Line {
	private Vertex p1;
	private Vertex p2;
	double xdiff, ydiff;
	double c, m, theta;

	public Line(Vertex p1, Vertex p2) {
		this.p1 = p1;
		this.p2 = p2;
		this.xdiff = p2.x - p1.x;
		this.ydiff = p2.y - p1.y;
		if (this.xdiff != 0) {
			this.m = this.ydiff / this.xdiff;
			this.c = ((p1.y + p2.y) - m * (p1.x + p2.x)) / 2;
			this.theta = Math.atan2(this.ydiff, this.xdiff);
		} else {
			if (this.ydiff > 0) {
				this.theta = Math.PI / 2.0;
			} else {
				this.theta = -Math.PI / 2.0;
			}
			this.c = 0;
		}
		
	}

	public String toString() {
		return "[" + this.p1.x + " " + this.p1.y + "; " + this.p2.x + " "
				+ this.p2.y + "]";
	}

	public double[] pIntersection(Vertex p, double eps, boolean precision) {

		double t[] = new double[2];
		/*
		 * Line 1: x1+(x2-x1)*t = x; x1+xdiff*t = x y1+(y2-y1)*t = y; y1+ydiff*t
		 * = y
		 * 
		 * (x-a)^2+(y-b)^2=eps^2 (x1+xdiff*t-a)^2+(y1+ydiff*t-b)^2=eps^2
		 * (xdiff^2+ydiff^2)t^2 + (2(x1-a)xdiff+2(y1-b)ydiff)t +
		 * (x1-a)^2+(y1-b)^2-eps^2=0
		 * 
		 * t = (-(2(x1-a)xdiff+2(y1-b)ydiff)) +-
		 * sqrt((2(x1-a)xdiff+2(y1-b)ydiff))^2 -
		 * 4(xdiff^2+ydiff^2)(x1^2+y1^2-eps^2)))/(2*(xdiff^2+ydiff^2))
		 */

		double b_t = 2 * ((p1.x - p.x) * this.xdiff + (p1.y - p.y) * this.ydiff);
		double c_t = (p1.x - p.x) * (p1.x - p.x) + (p1.y - p.y) * (p1.y - p.y)
				- eps * eps;
		double a_t = this.xdiff * this.xdiff + this.ydiff * this.ydiff;

		double newEps = eps + 0.01;
		double new_c_t = (p1.x - p.x) * (p1.x - p.x) + (p1.y - p.y)
				* (p1.y - p.y) - newEps * newEps;

		if (a_t == 0) {
			/*
			 * System.out.println("a_t="+a_t+",b_t="+b_t+", c_t="+c_t);
			 * System.out.println("disc="+ (b_t*b_t - 4*a_t*c_t));
			 * System.out.println(p1.toString()+p2.toString());
			 */
			return null;
		}

		double determinant = b_t * b_t - 4 * a_t * c_t;
		double newDeterminant = b_t * b_t - 4 * a_t * new_c_t;

		if (determinant >= 0) {
			t[1] = (-b_t + Math.sqrt(determinant)) / (2 * a_t);
			t[0] = (-b_t - Math.sqrt(determinant)) / (2 * a_t);
		} else if (precision && newDeterminant >= 0) {
			t[1] = (-b_t + Math.sqrt(0)) / (2 * a_t);
			t[0] = (-b_t - Math.sqrt(0)) / (2 * a_t);
		} else {
			return null;
		}
		return t;
	}

	public boolean pIntersection(Edge e, double eps) {

		/*
		 * Line 1: x1+(x2-x1)*t = x y1+(y2-y1)*t = y Line 2: y = mx + c
		 * y1+(y2-y1)*t = (x1+(x2-x1)*t)*m + c (y2-y1)*t - (x2-x1)*t*m = x1*m +
		 * c- y1 t = (x1*m + c - y1)/((y2-y1)-(x2-x1)*m)
		 */

		Line vline = e.line;

		Line line[] = this.getEpsilonNeiborhood(vline, eps);
		double startIntervalOnCurve, endIntervalOnCurve, startIntervalOnEdge, endIntervalOnEdge;
		double vstart = -1, vend = -1;

		if (Math.abs(this.theta) == Math.abs(line[0].theta)) {
			double t[] = this.getTParallel(vline, eps);
			if (t == null) {
				return false;
			}
			startIntervalOnCurve = t[0];
			endIntervalOnCurve = t[1];
			startIntervalOnEdge = 0;
			endIntervalOnEdge = 1;

		} else {

			if (line[0].xdiff == 0) {
				startIntervalOnCurve = this.m * line[0].p1.x + this.c;
				endIntervalOnCurve = this.m * line[1].p1.x + this.c;

				startIntervalOnEdge = (startIntervalOnCurve - line[0].p1.y)
						/ (line[0].ydiff);
				endIntervalOnEdge = (endIntervalOnCurve - line[1].p1.y)
						/ (line[1].ydiff);

				startIntervalOnCurve = (startIntervalOnCurve - this.p1.y)
						/ (this.ydiff);
				endIntervalOnCurve = (endIntervalOnCurve - this.p1.y)
						/ (this.ydiff);
			} else if (this.xdiff == 0) {
				startIntervalOnCurve = line[0].m * this.p1.x + line[0].c;
				endIntervalOnCurve = line[1].m * this.p1.x + line[1].c;

				startIntervalOnEdge = (startIntervalOnCurve - line[0].p1.y)
						/ (line[0].ydiff);
				endIntervalOnEdge = (endIntervalOnCurve - line[1].p1.y)
						/ (line[1].ydiff);

				startIntervalOnCurve = (startIntervalOnCurve - this.p1.y)
						/ (this.ydiff);
				endIntervalOnCurve = (endIntervalOnCurve - this.p1.y)
						/ (this.ydiff);
			} else {
				startIntervalOnCurve = -(this.c - line[0].c)
						/ (this.m - line[0].m);
				endIntervalOnCurve = -(this.c - line[1].c)
						/ (this.m - line[1].m);

				startIntervalOnEdge = (startIntervalOnCurve - line[0].p1.x)
						/ (line[0].p2.x - line[0].p1.x);
				endIntervalOnEdge = (endIntervalOnCurve - line[1].p1.x)
						/ (line[1].p2.x - line[1].p1.x);

				startIntervalOnCurve = (startIntervalOnCurve - this.p1.x)
						/ (this.p2.x - this.p1.x);
				endIntervalOnCurve = (endIntervalOnCurve - this.p1.x)
						/ (this.p2.x - this.p1.x);
			}
			/*if (line[0].p2.x - line[0].p1.x == 0) {

				return false;
			}*/

		}

		double intersectionWithDiscOne[] = new double[2];
		double intersectionWithDiscTwo[] = new double[2];
		double interval[] = new double[2];

		if (startIntervalOnCurve > endIntervalOnCurve) {
			double temp = startIntervalOnCurve;
			startIntervalOnCurve = endIntervalOnCurve;
			endIntervalOnCurve = temp;

			temp = startIntervalOnEdge;
			startIntervalOnEdge = endIntervalOnEdge;
			endIntervalOnEdge = temp;
		}

		intersectionWithDiscOne = this.pIntersection(vline.p1, eps, false);

		intersectionWithDiscTwo = this.pIntersection(vline.p2, eps, false);

		double min1 = 0, min2 = 0, max1 = 1, max2 = 1;

		// case one

		if (intersectionWithDiscOne == null && intersectionWithDiscTwo == null) {
			// case one
			if (startIntervalOnCurve > 1 || endIntervalOnCurve < 0)
				return false;
			// case two
			if ((startIntervalOnEdge > 1 && endIntervalOnEdge > 1)
					|| (startIntervalOnEdge < 0 && endIntervalOnEdge < 0))
				return false;

			if (startIntervalOnEdge >= 0 && startIntervalOnEdge <= 1
					&& endIntervalOnEdge >= 0 && endIntervalOnEdge <= 1) {
				startIntervalOnCurve = Math.max(0, startIntervalOnCurve);
				endIntervalOnCurve = Math.min(1, endIntervalOnCurve);
			} else
				return false;
		}

		// case two

		if (intersectionWithDiscOne != null) {

			min1 = Math.min(intersectionWithDiscOne[0], intersectionWithDiscOne[1]);
			max1 = Math.max(intersectionWithDiscOne[0], intersectionWithDiscOne[1]);

			if (((startIntervalOnEdge > 1 && endIntervalOnEdge > 1) || (startIntervalOnEdge < 0 && endIntervalOnEdge < 0))
					&& (min1 > 1 || max1 < 0))
				return false;

			if (startIntervalOnEdge < 0) {
				if (min1 <= 1)
					startIntervalOnCurve = Math.max(startIntervalOnCurve, min1);
				if (startIntervalOnCurve == min1)
					vstart = 0;
			}

			if (endIntervalOnEdge < 0) {
				if (max1 >= 0)
					endIntervalOnCurve = Math.min(endIntervalOnCurve, max1);
				if (endIntervalOnCurve == max1)
					vend = 0;
			}

			

		}
		// case three
		if (intersectionWithDiscTwo != null) {
			min2 = Math.min(intersectionWithDiscTwo[0], intersectionWithDiscTwo[1]);
			max2 = Math.max(intersectionWithDiscTwo[0], intersectionWithDiscTwo[1]);

			if (((startIntervalOnEdge > 1 && endIntervalOnEdge > 1) || (startIntervalOnEdge < 0 && endIntervalOnEdge < 0))
					&& (min2 > 1 || max2 < 0))
				return false;

			if (startIntervalOnEdge > 1) {
				if (min2 <= 1)
					startIntervalOnCurve = Math.max(startIntervalOnCurve, min2);
				if (startIntervalOnCurve == min2)
					vstart = 1;
			}

			if (endIntervalOnEdge > 1) {
				if (max2 >= 0)
					endIntervalOnCurve = Math.min(endIntervalOnCurve, max2);
				if (endIntervalOnCurve == max2)
					vend = 1;
			}
			
		}

		if (intersectionWithDiscOne == null && intersectionWithDiscTwo != null) {
			if (startIntervalOnEdge >= 0 && startIntervalOnEdge <= 1
					&& endIntervalOnEdge >= 0 && endIntervalOnEdge <= 1) {

			} else {
				if (startIntervalOnEdge >= 0 && startIntervalOnEdge <= 1) {
					if (startIntervalOnCurve > 1 && min2 > 1)
						return false;
					if (startIntervalOnCurve < 0 && max2 < 0)
						return false;
				}
				if (endIntervalOnEdge >= 0 && endIntervalOnEdge <= 1) {
					if (endIntervalOnCurve > 1 && min2 > 1)
						return false;
					if (endIntervalOnCurve < 0 && max2 < 0)
						return false;
				}
			}
		}
		if (intersectionWithDiscOne != null && intersectionWithDiscTwo == null) {
			if (startIntervalOnEdge >= 0 && startIntervalOnEdge <= 1
					&& endIntervalOnEdge >= 0 && endIntervalOnEdge <= 1) {

			} else {

				if (startIntervalOnEdge >= 0 && startIntervalOnEdge <= 1) {
					if (startIntervalOnCurve > 1 && min1 > 1)
						return false;
					if (startIntervalOnCurve < 0 && max1 < 0)
						return false;
				}
				if (endIntervalOnEdge >= 0 && endIntervalOnEdge <= 1) {
					if (endIntervalOnCurve > 1 && min1 > 1)
						return false;
					if (endIntervalOnCurve < 0 && max1 < 0)
						return false;
				}
			}
		}

		if (startIntervalOnCurve > endIntervalOnCurve) {
			double temp = startIntervalOnCurve;
			startIntervalOnCurve = endIntervalOnCurve;
			endIntervalOnCurve = temp;

			temp = vend;
			vend = vstart;
			vstart = temp;
		}
		// case one
		if (startIntervalOnCurve > 1 || endIntervalOnCurve < 0)
			return false;
		// case two
		
		if (Math.max(max1, max2) < 0 || Math.min(min1, min2) > 1)
			return false;

		double in1[] = new double[2];
		double in2[] = new double[2];
		try {
			interval[0] = Math.max(0, startIntervalOnCurve);

			if (vstart == -1) {

				in1 = vline.pIntersection(this.getVertex(interval[0]), eps,
						true);
				
				if (in1[0] >= 0 && in1[0] <= 1 && in1[1] >= 0 && in1[1] <= 1)
					vstart = (in1[0] + in1[1]) / 2;
				else if (in1[0] >= 0 && in1[0] <= 1)
					vstart = in1[0];
				else if (in1[1] >= 0 && in1[1] <= 1)
					vstart = in1[1];
				else
					vstart = 0;
			}

			interval[1] = Math.min(1, endIntervalOnCurve);

			if (vend == -1) {

				in2 = vline.pIntersection(this.getVertex(interval[1]), eps,
						true);
				
				if (in2[0] >= 0 && in2[0] <= 1 && in2[1] >= 0 && in2[1] <= 1)
					vend = (in2[0] + in2[1]) / 2;
				else if (in2[0] >= 0 && in2[0] <= 1)
					vend = in2[0];
				else if (in2[1] >= 0 && in2[1] <= 1)
					vend = in2[1];
				else
					vend = 1;
			}
			e.cstart = interval[0];
			e.cend = interval[1];
			e.vstart = vstart;
			e.vend = vend;

			return true;
		} catch (Exception ex) {
			System.out.println("eline: " + e.line.toString() + " " + e.line.m
					+ " " + e.line.theta);
			System.out.println("this: " + this.toString() + " " + this.m + " "
					+ this.theta);
			System.out.println(startIntervalOnCurve + " " + endIntervalOnCurve
					+ " " + startIntervalOnEdge + " " + endIntervalOnEdge);
			System.out.println(min1 + " " + max1 + " " + min2 + " " + max2);
			System.out.println(this.getVertex(interval[0]).toString());
			System.out.println(this.getVertex(interval[1]).toString());
			
			System.out.println(ex.toString());
			System.exit(0);
			return false;
		}
	}

	public Vertex getVertex(double t) {
		return new Vertex(this.p1.x + this.xdiff * t, this.p1.y + this.ydiff
				* t);
	}

	public Line[] getEpsilonNeiborhood(Line vline, double eps) {

		Line line[] = new Line[2];
		double dx, dy;
		dx = eps * Math.cos(vline.theta + Math.PI / 2);
		dy = eps * Math.sin(vline.theta + Math.PI / 2);
		// dx_l = eps*Math.cos(vline.theta + 3*Math.PI/2);
		// dy_l = eps*Math.sin(vline.theta + 3*Math.PI/2);

		line[0] = new Line(new Vertex(vline.p1.x - dx, vline.p1.y - dy),
				new Vertex(vline.p2.x - dx, vline.p2.y - dy));
		line[1] = new Line(new Vertex(vline.p1.x + dx, vline.p1.y + dy),
				new Vertex(vline.p2.x + dx, vline.p2.y + dy));

		if (line[0].theta != line[1].theta) {

			line[0].m = line[1].m;
			line[0].theta = line[1].theta;
		}
		return line;
	}

	public double[] getTParallel(Line line, double eps) {
		double t[] = new double[2];
		double newm;
		double x1, y1, x2, y2;
		if (Math.abs(line.theta) == Math.PI / 2) {
			newm = 0;
			x1 = this.p1.x;
			y1 = line.p1.y;
			x2 = this.p2.x;
			y2 = line.p2.y;
		} else if (Math.abs(line.theta) == 0) {
			newm = 0;
			x1 = line.p1.x;
			y1 = this.p1.y;
			x2 = line.p2.x;
			y2 = this.p2.y;
		} else {
			newm = 1 / line.m;
			double c1 = line.p1.y + newm * line.p1.x;
			double c2 = line.p2.y + newm * line.p2.x;

			x1 = (c1 - this.c) / (this.m + newm);
			y1 = this.m * x1 + this.c;

			x2 = (c2 - this.c) / (this.m + newm);
			y2 = this.m * x2 + this.c;
		}

		if (Math.sqrt((line.p1.x - x1) * (line.p1.x - x1) + (line.p1.y - y1)
				* (line.p1.y - y1)) > eps)
			return null;

		if (this.xdiff != 0) {
			t[0] = (x1 - this.p1.x) / this.xdiff;
			t[1] = (x2 - this.p1.x) / this.xdiff;
		} else {
			t[0] = (y1 - this.p1.y) / this.ydiff;
			t[1] = (y2 - this.p1.y) / this.ydiff;
		}

		double min = Math.min(t[0], t[1]);
		double max = Math.max(t[0], t[1]);

		if (min > 1 || max < 0)
			return null;
		t[0] = min;
		t[1] = max;

		return t;
	}

	public double distance(Vertex p) {

		double b_t = 2 * ((p1.x - p.x) * this.xdiff + (p1.y - p.y) * this.ydiff);
		double c_t = (p1.x - p.x) * (p1.x - p.x) + (p1.y - p.y) * (p1.y - p.y);// -
																				// eps*eps;
		double a_t = this.xdiff * this.xdiff + this.ydiff * this.ydiff;
		double disc = Math.sqrt((b_t * b_t - 4 * a_t * c_t) / (4 * a_t));
		double t[] = this.pIntersection(p, disc, false);
		if (t == null || t[0] > 1 || t[0] < 0)
			return Math.min(
					Math.sqrt((p1.x - p.y) * (p1.x - p.y) + (p1.y - p.y)
							* (p1.y - p.y)),
					Math.sqrt((p2.x - p.y) * (p2.x - p.y) + (p2.y - p.y)
							* (p2.y - p.y)));
		else
			return disc;
	}

	
}
