#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <vector>
#include <map>
#include <utility>
#include <math.h>

using namespace std;

int WINDOWSIZE = 750;

/* isSimple
 * 
 * Inputs:
 *     A vector of points which form a polygon.
 *
 * Return Values:
 *     This function returns a bool indicating whether or not the inputted polygon is simple.
 *
 * Functional Description:
 *     This function first creates all of the polygon edges from the inputted polygon points. This function then executes
 *     a double for loop through the polygon edges to see if any of the edges intersect. If no edges intersect, this function
 *     returns true; otherwise, it returns false.
 */

bool isSimple(vector<point2D> polygonPoints) {
  
	vector<lineSegment2D> polygonEdges = makePolygonEdges(polygonPoints);

	for (int i = 0; i < polygonEdges.size(); i++) {
		for (int j = 0; j < polygonEdges.size(); j++) {
			if (properIntersect(polygonEdges[i].p1, polygonEdges[i].p2, polygonEdges[j].p1, polygonEdges[j].p2)) {
		 		return false;
		 	}
		}
	} 

	return true;
}

/* polygonShift
 * 
 * Inputs:
 *     A vector of points which form a polygon.
 *
 * Return Values:
 *     A vector of shifted polygon points, where the guardPoint is on the origin.
 *
 * Functional Description:
 *     This function iterates through all the points in the polygon and shifts all the vertices by both the x and y coordinates
 *     of the origin, resulting in the guardPoint being the origin.
 */

vector<point2D> polygonShift(vector<point2D> originalPolygonPoints, point2D guardPoint) {

	for (int i = 0; i < originalPolygonPoints.size(); i++) {
		originalPolygonPoints[i].x -= guardPoint.x;
		originalPolygonPoints[i].y -= guardPoint.y;
	} 

	return originalPolygonPoints;

}

/* computeXIntersection
 * 
 * Inputs:
 *     The x and y coordinates of two points which form a line.
 *
 * Return Values:
 *     The x-coordinate intersection between the line which is an argument and the x-axis.
 *
 */

double computeXIntersection (int firstPointX, int firstPointY, int secondPointX, int secondPointY) {
	return ((firstPointX * secondPointY - secondPointX * firstPointY)) / ((double)(secondPointY - firstPointY));
}

/* isInPolygon
 * 
 * Inputs:
 *     - A vector points representing the polygon in question.
 * 	   - A point representing the current position of the guard in the polygon.
 *
 * Return Values:
 *     This function returns a bool indicating if the guardPoint is in the polygon argument or not.
 *
 * Functional Description:
 *     
 */

bool isInPolygon(vector<point2D> polygonPoints, point2D guardPoint) { //Need last mouse click for point

	int previousPoint; // = i1
	int numberOfPoints = polygonPoints.size();
	double xAxisIntersection;
	int Rcross = 0;
	int Lcross = 0;
	bool Rstrad, Lstrad;

	vector<point2D> shiftedPoints;

	shiftedPoints = polygonShift(polygonPoints, guardPoint); //Points are shifted by guardPoint

	for (int i = 0; i < polygonPoints.size(); i++) { // i = first edge point

		if (shiftedPoints[i].x == 0 && shiftedPoints[i].y == 0) {
			return true;
		}

		previousPoint = (i + numberOfPoints - 1) % numberOfPoints;

		Rstrad = (shiftedPoints[i].y > 0) != (shiftedPoints[previousPoint].y > 0);
		Lstrad = (shiftedPoints[i].y < 0) != (shiftedPoints[previousPoint].y < 0);

		if (Rstrad || Lstrad) {
			xAxisIntersection = computeXIntersection(shiftedPoints[i].x, shiftedPoints[i].y, shiftedPoints[previousPoint].x, shiftedPoints[previousPoint].y);

			if (Rstrad && xAxisIntersection > 0) {
				Rcross += 1;
			}

			if (Lstrad && xAxisIntersection > 0) {
				Lcross += 1;
			}

		}

	}

	if ((Rcross % 2) != (Lcross % 2)) {
		return true;
	}

	if ( (Rcross % 2) == 1) {
		return true;
	}

	else {
		return false;
	}

}

/* collinear
 * 
 * Inputs:
 *     Three point structs.
 *
 * Return Values:
 *     A bool indicating whether or not the three inputted points are collinear or not.
 *
 */

bool collinear(point2D a, point2D b, point2D c) {
	return (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) == 0;
}

/* Area2
 * 
 * Inputs:
 *     Three point structs.
 *
 * Return Values:
 *     An int representing the area formed by the three points as arguments.
 *
 */

int Area2(point2D a, point2D b, point2D c) {
	return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}

/* Xor
 * 
 * Inputs:
 *     Two generic bools.
 *
 * Return Values:
 *     A bool indicating if the two inputted bools differ.
 *
 */

bool Xor(bool x, bool y) {
	return !x ^ !y;
}

/* left
 * 
 * Inputs:
 *     Three points structs.
 *
 * Return Values:
 *     A bool indicating whether of not point c is to the left of points a and b.
 *
 */

bool left(point2D a, point2D b, point2D c) {
	return Area2(a,b,c) > 0;
}

/* properIntersect
 * 
 * Inputs:
 *     Four points which form two individual lines.
 *
 * Return Values:
 *     A bool indicating whether or not the lines have a proper intersection.
 *
 */

bool properIntersect(point2D a, point2D b, point2D c, point2D d) {
	if (collinear(a,b,c) || collinear(a,b,d) || collinear(c,d,a) || collinear(c,d,b)) {
		return false;
	}

	return Xor(left(a,b,c), left(a,b,d)) && Xor(left(c,d,a), left(c,d,b));
}

/* between
 * 
 * Inputs:
 *     Three points structs.
 *
 * Return Values:
 *     A bool indicating whether or not point c is between points a and b.
 *
 */

bool between(point2D a, point2D b, point2D c) {
	
	if (!collinear(a,b,c)) {
		return false;
	}

	if (a.x != b.x) {
		return ((a.x <= c.x) && (c.x <= b.x)) || ((a.x >= c.x) && (c.x >= b.x));
	}

	else {
		return ((a.y <= c.y) && (c.y <= b.y)) || ((a.y >= c.y) && (c.y >= b.y));
	}

}

/* improperIntersect
 * 
 * Inputs:
 *     Four points which form two individual lines.
 *
 * Return Values:
 *     A bool indicating whether or not the lines have an improper intersection.
 *
 */

bool improperIntersect(point2D a, point2D b, point2D c, point2D d) {
	return between(a, b, c) || between(a, b, d) || between(c, d, a) || between(c, d, b);
}

/* segmentLength
 * 
 * Inputs:
 *     Two points which form a line segment.
 *
 * Return Values:
 *     The length of the segment between the two points.
 *
 */

double segmentLength(point2D point1, point2D point2) {
	return sqrt((pow((point1.x - point2.x), 2)) + (pow((point1.y - point2.y), 2)));
}

/* extendLine
 * 
 * Inputs:
 *     - A point which is one of the verticies in a polygon.
 *     - A point which represents a guard.
 *
 * Return Values:
 *     This function returns a line segment from the guard point passing through the vertex to the edge of the window.
 *
 * Functional Description:
 *     This function generates a new end point for the segment by using the current slope between the guard point and vertex,
 *     and extending this slope to the edge of the window.
 *     
 */

lineSegment2D extendLine(point2D polygonPoint, point2D guardPoint) {

	lineSegment2D extendedLine;
	point2D newEndPoint;

	newEndPoint.x = WINDOWSIZE*(polygonPoint.x - guardPoint.x) + polygonPoint.x;
	newEndPoint.y = WINDOWSIZE*(polygonPoint.y - guardPoint.y) + polygonPoint.y;

	extendedLine.p1 = guardPoint;
	extendedLine.p2 = newEndPoint;

	return extendedLine;

}

/* extendLineBothDirections
 * 
 * Inputs:
 *     - A point which is one of the verticies in a polygon.
 *     - A point which represents a guard.
 *
 * Return Values:
 *     This function returns a line segment passing through both the guard point and vertex, extending to both ends of the window.
 *
 * Functional Description:
 *     This function generates a new start and end point for the segment by using the current slope between the guard point and vertex,
 *     and extending this slope to both the edges of the window.
 *     
 */

lineSegment2D extendLineBothDirections(point2D polygonPoint, point2D guardPoint) {

	lineSegment2D extendedLine;
	point2D newEndPoint, newEndPoint2;

	newEndPoint.x = WINDOWSIZE*(polygonPoint.x - guardPoint.x) + polygonPoint.x;
	newEndPoint.y = WINDOWSIZE*(polygonPoint.y - guardPoint.y) + polygonPoint.y;

	newEndPoint2.x = WINDOWSIZE*(guardPoint.x - polygonPoint.x) + guardPoint.x;
	newEndPoint2.y = WINDOWSIZE*(guardPoint.y - polygonPoint.y) + guardPoint.y;

	extendedLine.p1 = newEndPoint;
	extendedLine.p2 = newEndPoint2;

	return extendedLine;

}

/* linesIntersect
 * 
 * Inputs:
 *     - The extended line between the guard point and a point on the polygon.
 *
 * Return Values:
 *     A bool indicating whether or not the two lines intersect properly or improperly.
 *
 */

bool linesIntersect(lineSegment2D guardVertexEdge, lineSegment2D polygonEdge) {

	if (properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2) || improperIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2)) {
		return true;
	}

	return false;

}

/* grazePoint
 * 
 * Inputs:
 *     - The point before a vertex in the polygon.
 *     - The point after a vertex in the polygon.
 * 	   - The extended line between a guard point and the vertex in the polygon.
 *
 * Return Values:
 *     A bool indicating whether or not a vertex in a polygon is a graze point, i.e. the guard can see past this point.
 *
 * Functional Description:
 *     This function determines if there is a proper intersection between a segment between the points before and after a point
 *     in the polygon, and the extended line between a guard point and the vertex in the polygon. If there is a proper intersection between
 *     these lines, then that vertex is not a graze point. Otherwise, it is.
 *     
 */

bool grazePoint(point2D polygonPointV1, point2D polygonPointV2, lineSegment2D guardVertexEdge) {

	lineSegment2D prevNextLine;

	prevNextLine.p1 = polygonPointV1;
	prevNextLine.p2 = polygonPointV2;

	if (properIntersect(prevNextLine.p1, prevNextLine.p2, guardVertexEdge.p1, guardVertexEdge.p2)) {
		return false;
	}

	else {
		return true;
	}

}

/* intersectionPoint
 * 
 * Inputs:
 *     Four point structs which form two individual lines.
 *
 * Return Values:
 *     The point of intersection between the two lines.
 *
 */

point2D intersectionPoint(point2D a, point2D b, point2D c, point2D d) {

	double s;
	double num, denom;
	point2D intersectionPoint;

	denom = a.x * (double) (d.y - c.y) + b.x * (double) (c.y - d.y) + d.x * (double) (b.y - a.y) + c.x * (double) (a.y - b.y);

	num = a.x * (double) (d.y - c.y) + c.x * (double) (a.y - d.y) + d.x * (double) (c.y - a.y);

	s = num / denom;

	intersectionPoint.x = a.x + s * (b.x - a.x);
	intersectionPoint.y = a.y + s * (b.y - a.y);

	return intersectionPoint;

}

/* eventSorterDistance
 * 
 * Inputs:
 *     Two event structs.
 *
 * Return Values:
 *     A bool indicating which event's distance from the guard is less. This is used to sort events based on the event's
 * 	   distance from the guard.
 *
 */

bool eventSorterDistance(const eventPoint& pointOne, const eventPoint& pointTwo) {
	return pointOne.distanceFromGuard < pointTwo.distanceFromGuard;
}

/* pointDistance
 * 
 * Inputs:
 *     Two individual points.
 *
 * Return Values:
 *     The distance between two points.
 *
 */

double pointDistance(point2D event, point2D guardPoint) {
	return sqrt((pow((event.x - guardPoint.x), 2)) + (pow((event.y - guardPoint.y), 2)));
}

/* pointExists
 * 
 * Inputs:
 *     - A vector containing the current visible polygon points.
 *     - A new point to be inserted into the visible polygon points.
 *
 * Return Values:
 *     This function returns a bool indicating whether or not the point to be inserted already exists, in which case,
 * 	   this point should not be inserted into the vector again.
 *
 * Functional Description:
 *     This function iterates through the current vector of visible points, and checks every point to see if the point
 *     to be inserted already exists.
 *     
 */

bool pointExists(vector<point2D> visiblePolygon, point2D newPoint) {

	for (int i = 0; i < visiblePolygon.size(); i++) {
		if (visiblePolygon[i].x == newPoint.x && visiblePolygon[i].y == newPoint.y) {
			return true;
		}
	}

	return false;

}

/* eventExists
 * 
 * Inputs:
 *     - A vector containing the current events.
 *     - A new event to be inserted into the events vector.
 *
 * Return Values:
 *     This function returns a bool indicating whether or not the event to be inserted already exists, in which case,
 * 	   this event should not be inserted into the vector again.
 *
 * Functional Description:
 *     This function iterates through the current vector of events, and checks every event to see if the event to be inserted
 *     already exists.
 *     
 */

bool eventExists(vector<eventPoint> events, point2D newPoint) {

	for (int i = 0; i < events.size(); i++) {
		if (events[i].point.x == newPoint.x && events[i].point.y == newPoint.y) {
			return true;
		}
	}

	return false;

}

/* makePolygonEdges
 * 
 * Inputs:
 *     A vector of all points in a polygon. The points are sorted so that the neighbouring points in the vector are
 *     neighbouring points in the polygon.
 *
 * Return Values:
 *     This function returns the vector of line segments between neighbouring points in the polygon.
 *
 * Functional Description:
 *     This function iterates through all the points in polygon, and creates a new line segment between the current point, and the
 *     next point in the vector. This results in line segments between all vertices in the polygon.
 */

vector<lineSegment2D> makePolygonEdges(vector<point2D> polygonPoints) {

	vector<lineSegment2D> polygonEdges;
	lineSegment2D newEdge;

	for (int i = 0; i < polygonPoints.size() - 1; i++) {
		newEdge.p1 = polygonPoints[i];
		newEdge.p2 = polygonPoints[i+1];
		polygonEdges.push_back(newEdge);
	}

	newEdge.p1 = polygonPoints.back();
	newEdge.p2 = polygonPoints.front();
	polygonEdges.push_back(newEdge);

	return polygonEdges;

}

/* pointComparer
 * 
 * Inputs:
 *     Two point structs
 *
 * Return Values:
 *     A bool indicating whether or not the two points are the same point.
 *
 */

bool pointComparer(const point2D& pointOne, const point2D& pointTwo) {

  if (pointOne.x == pointTwo.x && pointOne.y == pointTwo.y) {
    return true;
  }

  return false;

}

/* distanceSorter
 * 
 * Inputs:
 *     Two point structs
 *
 * Return Values:
 *     A bool indicating which point's distance from the guard is less. This is used to sort points based on the point's
 * 	   distance from the guard.
 *
 */

bool distanceSort(const point2D& pointOne, const point2D& pointTwo) {
  return pointOne.distanceFromPreviousPoint < pointTwo.distanceFromPreviousPoint;
}

/* remove_coincident_points
 * 
 * Inputs:
 *     A pointer to a vector of points.
 *
 * Return Values:
 *     N/A
 *
 * Functional Description:
 *     This function removes all conincident points from the vector being pointed to.
 */

void remove_coincident_points(vector<point2D> *points) {
  vector<point2D>::iterator unique_end = (unique((*points).begin(), (*points).end(), pointComparer));
  (*points).erase(unique_end, (*points).end());
  
}

/* reformPolygon
 * 
 * Inputs:
 *     - A vector of visible polygon points with respect to the guard point.
 * 	   - A vector of line segments representing the edges of the polygon.
 *
 * Return Values:
 *     This function returns a vector of points representing the visible polygon with respect to the guard point. This vector
 *     is sorted based on the order in which the visible points appear along the edges of the polygon.
 *
 * Functional Description:
 *     This function iterates through the ordered polygon edges, and then checks all the visible points that appear on that edge (
 *     these points are pushed onto a vector). If there are no points on that edge, the function continues onto the next edge. If there
 *     is only one point on that edge, that point is added to the ordered vector of points. If there is more than one point on that edge,
 *     then the points are sorted based on their distance relative to the first point of the polygon edge; these points are then inserted
 *     into the sorted array. After all edges have been processed, the resulting vector is the vector of ordered points.
 */

vector<point2D> reformPolygon(vector<point2D> visiblePolygonPoints, vector<lineSegment2D> polygonEdges) {

	vector<point2D> edgePoints, reformedPolygon;

	for (int i = 0; i < polygonEdges.size(); i++) {  //Ordered polygon edges.
		for (int j = 0; j < visiblePolygonPoints.size(); j++) {
			if (visiblePolygonPoints[j].intersectionEdge == i) {
				edgePoints.push_back(visiblePolygonPoints[j]);
			}
		}

		if (edgePoints.size() == 0) { //No visible points on that edge.
			continue;
		}

		else if (edgePoints.size() == 1) { //Only one visible point on that edge.
			reformedPolygon.push_back(edgePoints[0]);
		}

		else { 
			sort(edgePoints.begin(), edgePoints.end(), distanceSort); 

			for (int j = 0; j < edgePoints.size(); j++) {
				reformedPolygon.push_back(edgePoints[j]);
			}
		}

		edgePoints.clear();

	}

	return reformedPolygon;

}

/* initializeNewEvent
 * 
 * Inputs:
 *     - A pointer to a new event
 *     - An int indicating which edge the event appears on.
 *     - A lineSegment which the points appears on.
 *     - The extended line segment between the guard and a point on the polygon
 *     - A point representing the guard point.
 *
 * Return Values:
 *     N/A
 *
 * Functional Description:
 *     This function initializes all the properties of a new event struct.
 */

void initializeNewEvent(eventPoint *newEvent, int intersectionEdge, lineSegment2D polygonEdge, lineSegment2D guardVertexEdge, point2D guardPoint) {
	(*newEvent).intersectionEdge = intersectionEdge;
	(*newEvent).distanceFromPreviousPoint = pointDistance((*newEvent).point, polygonEdge.p1);
	(*newEvent).intersectionType = properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2);
	(*newEvent).distanceFromGuard = pointDistance((*newEvent).point, guardPoint);
}

/* initializeNewPoint
 * 
 * Inputs:
 *     - A pointer to a new event
 *     - A point representing the guard point.
 *     - A point which is the point before the point being initializes on the polygon.
 *
 * Return Values:
 *     N/A
 *
 * Functional Description:
 *     This function initializes the point component of the event struct.
 */

void initializePoint(eventPoint *newEvent, point2D guardPoint, point2D prevPolygonPoint) {
	(*newEvent).point.distanceFromGuard = (*newEvent).distanceFromGuard;
	(*newEvent).point.intersectionType = (*newEvent).intersectionType;
	(*newEvent).point.left = left(guardPoint, (*newEvent).point, prevPolygonPoint);
	(*newEvent).point.intersectionEdge = (*newEvent).intersectionEdge;
	(*newEvent).point.distanceFromPreviousPoint = (*newEvent).distanceFromPreviousPoint;
}

/* visiblePolygon
 * 
 * Inputs:
 *     - A vector points representing the polygon in question.
 * 	   - A point representing the current position of the guard in the polygon.
 *
 * Return Values:
 *     This function returns a vector of all the points on the polygon which the guard can see.
 *
 * Functional Description:
 *     This function first removes all coincident points in the polygon vector, before generating all of the edges in this polygon.
 *     This function then iterates through all of the points in the polygon, and creates the extended line between the guard point,
 *     and the current point in the polygon. Then this function determines all the intersection points between that extended line,
 *     and all edges in the polygon. If any intersection takes place, a new corresponding event is created. These events are then
 *     sorted based on their distance from the guard point. This function then iterates through the sorted events, and checks what
 *     kind of intersection the events are. If the event is an improper intersection and a graze point, this point is added to the
 *     vector of visible polygon points, and the function continues through the vector of events (as we can see past the graze point,
 *     and may be able to see more points along the same line). Otherwise, if the event is a proper intersection or it is not a graze
 *     point, the event is pushed back, before the function breaks the loop through the events and continues onto the next polygon vertex
 *     (as we cannot see past that point). This function finally reorders the visible polygon by calling the reformPolygon function.
 *     
 */

vector<point2D> visiblePolygon(vector<point2D> polygonPoints, point2D guardPoint) {

	remove_coincident_points(&polygonPoints);


	vector<point2D> visiblePolygon;
	vector<lineSegment2D> polygonEdges = makePolygonEdges(polygonPoints);
	vector<eventPoint> intersectionEvents, tempContainer;

	lineSegment2D guardVertexEdge;
	eventPoint newEvent;

	int prevPointIndex, nextPointIndex;

	for (int i = 0; i < polygonPoints.size(); i++) {

		if (i == 0) {
			prevPointIndex = polygonPoints.size() - 1;
			nextPointIndex = 1;
		}

		else if (i == polygonPoints.size() - 1) {
			prevPointIndex = polygonPoints.size() - 2;
			nextPointIndex = 0;
		}

		else {
			prevPointIndex = i - 1;
			nextPointIndex = i + 1;
		}

		guardVertexEdge = extendLine(polygonPoints[i], guardPoint);

		for (int j = 0; j < polygonEdges.size(); j++) {

			if (linesIntersect(guardVertexEdge, polygonEdges[j])) {
		 		newEvent.point = intersectionPoint(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdges[j].p1, polygonEdges[j].p2);

		 		if (eventExists(intersectionEvents, newEvent.point)) {
					continue;
				}

				initializeNewEvent(&newEvent, j, polygonEdges[j], guardVertexEdge, guardPoint);
		 		intersectionEvents.push_back(newEvent);

		 	}

		}

		sort(intersectionEvents.begin(), intersectionEvents.end(), eventSorterDistance);

		for (int k = 0; k < intersectionEvents.size(); k++) {

			initializePoint(&intersectionEvents[k], guardPoint, polygonPoints[prevPointIndex]);
			guardVertexEdge = extendLineBothDirections(polygonPoints[i], guardPoint);

			if (intersectionEvents[k].intersectionType == 0 && grazePoint(polygonPoints[prevPointIndex], polygonPoints[nextPointIndex], guardVertexEdge)) {
				if (!pointExists(visiblePolygon, intersectionEvents[k].point)) {
					visiblePolygon.push_back(intersectionEvents[k].point);
				}
			}

			else { 
				if (!pointExists(visiblePolygon, intersectionEvents[k].point)) {
					visiblePolygon.push_back(intersectionEvents[k].point);
				}

				break;
			}

		}

		intersectionEvents.clear();

	}

	visiblePolygon = reformPolygon(visiblePolygon, polygonEdges);

	return visiblePolygon;
}

/* triangulatedPolygon
 * 
 * Inputs:
 *     - A vector points representing the polygon in question.
 * 	   - A point representing the current position of the guard in the polygon.
 *
 * Return Values:
 *     This function returns a vector or vector of points of the individual triangles forming the visible polygon.
 *
 * Functional Description:
 *     This function iterates through all the points in the visible polygon vector, and creates a new vector of points
 *     (which represent an individual triangle), and pushes back the current point in the visible points vector, the next point
 *     in the vector, as well the guard point. As the guard point can see all vertices in the visible polygon vector, the guard
 *     point is used as a vertex for all triangles.
 */

vector<vector<point2D> > triangulatedPolygon(vector<point2D> polygonPoints, point2D guardPoint) {

	vector<vector<point2D> > triangulatedPolygon;
	vector<point2D> visiblePolygonPoints, singleTriangle;
	point2D p1, p2;

	visiblePolygonPoints = visiblePolygon(polygonPoints, guardPoint);

	for (int i = 0; i < visiblePolygonPoints.size() - 1; i ++) {
		p1 = visiblePolygonPoints[i];
		p2 = visiblePolygonPoints[i+1];

		singleTriangle.push_back(p1);
		singleTriangle.push_back(p2);
		singleTriangle.push_back(guardPoint);

		triangulatedPolygon.push_back(singleTriangle);

		singleTriangle.clear();

	}

	p1 = visiblePolygonPoints.back();
	p2 = visiblePolygonPoints.front();

	singleTriangle.push_back(p1);
	singleTriangle.push_back(guardPoint);
	singleTriangle.push_back(p2);

	triangulatedPolygon.push_back(singleTriangle);

	return triangulatedPolygon;

}