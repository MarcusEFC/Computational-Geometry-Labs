#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <vector>
#include <queue>
#include <map>
#include <utility>
#include <math.h>

using namespace std; 

int WINDOWSIZE = 750;

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

vector<point2D> polygonShift(vector<point2D> originalPolygonPoints, point2D guardPoint) {

	for (int i = 0; i < originalPolygonPoints.size(); i++) {
		originalPolygonPoints[i].x -= guardPoint.x;
		originalPolygonPoints[i].y -= guardPoint.y;
	} 

	return originalPolygonPoints;

}

double computeXIntersection (int firstPointX, int firstPointY, int secondPointX, int secondPointY) {
	return ((firstPointX * secondPointY - secondPointX * firstPointY)) / ((double)(secondPointY - firstPointY));

}

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
			return false;
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
		return false;
	}

	if ( (Rcross % 2) == 1) {
		return true;
	}

	else {
		return false;
	}

}

bool collinear(point2D a, point2D b, point2D c) {
	return (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) == 0;
}

int Area2(point2D a, point2D b, point2D c) {
	return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}

bool Xor(bool x, bool y) {
	return !x ^ !y;
}

bool left(point2D a, point2D b, point2D c) {
	return Area2(a,b,c) > 0;
}

bool properIntersect(point2D a, point2D b, point2D c, point2D d) {
	if (collinear(a,b,c) || collinear(a,b,d) || collinear(c,d,a) || collinear(c,d,b)) {
		return false;
	}

	return Xor(left(a,b,c), left(a,b,d)) && Xor(left(c,d,a), left(c,d,b));
}

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

bool improperIntersect(point2D a, point2D b, point2D c, point2D d) {
	return between(a, b, c) || between(a, b, d) || between(c, d, a) || between(c, d, b);
}

double segmentLength(point2D point1, point2D point2) {
	return sqrt((pow((point1.x - point2.x), 2)) + (pow((point1.y - point2.y), 2)));
}

lineSegment2D extendLine(point2D polygonPoint, point2D guardPoint) {

	lineSegment2D extendedLine;
	point2D newEndPoint;

	newEndPoint.x = WINDOWSIZE*(polygonPoint.x - guardPoint.x) + polygonPoint.x;
	newEndPoint.y = WINDOWSIZE*(polygonPoint.y - guardPoint.y) + polygonPoint.y;

	extendedLine.p1 = guardPoint;
	extendedLine.p2 = newEndPoint;

	return extendedLine;

}

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

bool linesIntersect(lineSegment2D guardVertexEdge, lineSegment2D polygonEdge) {

	if (properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2)) {
		return true;
	}

	else if (improperIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2)) {
		return true;
	}

	else {
		return false;
	}

}

double slopeOfLine(lineSegment2D line) {
	return (line.p1.y - line.p2.y) / (line.p1.x - line.p2.x);
}

double yIntercept(lineSegment2D line, double slope) {
	return line.p1.y - (slope * line.p1.x);
}

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

bool eventSorterDistance(const eventPoint& pointOne, const eventPoint& pointTwo) {
	return pointOne.distanceFromGuard < pointTwo.distanceFromGuard;
}

bool pointSorterDistance(const point2D& pointOne, const point2D& pointTwo) {
	return pointOne.distanceFromGuard < pointTwo.distanceFromGuard;
}

bool pointSorterDistanceStart(const point2D& pointOne, const point2D& pointTwo) {
	return pointOne.distanceFromStart < pointTwo.distanceFromStart;
}

bool pointSorterNegativeDistance(const point2D& pointOne, const point2D& pointTwo) {
	return pointOne.distanceFromGuard > pointTwo.distanceFromGuard;
}

double pointDistance(point2D event, point2D guardPoint) {
	return sqrt((pow((event.x - guardPoint.x), 2)) + (pow((event.y - guardPoint.y), 2)));
}

vector<eventPoint> returnEvents(vector<point2D> polygonPoints, point2D guardPoint) { //Returns polygon visible to guard.

	vector<point2D> visiblePolygon;
	vector<lineSegment2D> polygonEdges = makePolygonEdges(polygonPoints);
	vector<eventPoint> intersectionEvents;

	lineSegment2D guardVertexEdge;
	eventPoint newEvent;

	for (int i = 0; i < polygonPoints.size(); i++) { //Need to check if point is on a vertex. If so, just add it.

		guardVertexEdge = extendLine(polygonPoints[i], guardPoint);

		for (int j = 0; j < polygonEdges.size(); j++) {

			if (linesIntersect(guardVertexEdge, polygonEdges[j]) == true) {
		 		newEvent.point = intersectionPoint(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdges[j].p1, polygonEdges[j].p2);
		 		newEvent.intersectionType = properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdges[j].p1, polygonEdges[j].p2);
		 		newEvent.distanceFromGuard = pointDistance(newEvent.point, guardPoint);
		 		intersectionEvents.push_back(newEvent);

		 	}

		}

	}
	return intersectionEvents;

}

bool pointExists(vector<point2D> visiblePolygon, point2D newPoint) {

	for (int i = 0; i < visiblePolygon.size(); i++) {
		if (visiblePolygon[i].x == newPoint.x && visiblePolygon[i].y == newPoint.y) {
			return true;
		}
	}

	return false;

}

bool eventExists(vector<eventPoint> events, point2D newPoint) {

	for (int i = 0; i < events.size(); i++) {
		if (events[i].point.x == newPoint.x && events[i].point.y == newPoint.y) {
			return true;
		}
	}

	return false;

}

vector<lineSegment2D> makePolygonEdges(vector<point2D> polygonPoints) { //Returns all edges in polygon

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

double computedAngle(point2D referencePoint, point2D newPoint) {

	double computeAngle = atan2((newPoint.y - referencePoint.y), (newPoint.x - referencePoint.x));

	if (computeAngle < 0) {
        computeAngle = computeAngle + 2*M_PI;
    }

    return computeAngle;

}

bool pointComparer(const point2D& pointOne, const point2D& pointTwo) {

  if (pointOne.x == pointTwo.x && pointOne.y == pointTwo.y) {
    return true;
  }

  return false;

}

bool distanceSort(const point2D& pointOne, const point2D& pointTwo) {
  return pointOne.distanceFromPreviousPoint < pointTwo.distanceFromPreviousPoint;

}

void remove_coincident_points(vector<point2D> *points) {
  vector<point2D>::iterator unique_end = (unique((*points).begin(), (*points).end(), pointComparer));
  (*points).erase(unique_end, (*points).end());
  
}

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

void initializeNewEvent(eventPoint *newEvent, int intersectionEdge, lineSegment2D polygonEdge, lineSegment2D guardVertexEdge, point2D guardPoint) {
	(*newEvent).intersectionEdge = intersectionEdge;
	(*newEvent).distanceFromPreviousPoint = pointDistance((*newEvent).point, polygonEdge.p1);
	(*newEvent).intersectionType = properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2);
	(*newEvent).distanceFromGuard = pointDistance((*newEvent).point, guardPoint);

}

void initializePoint(eventPoint *newEvent, point2D guardPoint, point2D prevPolygonPoint) {
	(*newEvent).point.distanceFromGuard = (*newEvent).distanceFromGuard;
	(*newEvent).point.guardAngle = computedAngle(guardPoint, (*newEvent).point);
	(*newEvent).point.intersectionType = (*newEvent).intersectionType;
	(*newEvent).point.left = left(guardPoint, (*newEvent).point, prevPolygonPoint);
	(*newEvent).point.intersectionEdge = (*newEvent).intersectionEdge;
	(*newEvent).point.distanceFromPreviousPoint = (*newEvent).distanceFromPreviousPoint;
}

vector<point2D> visiblePolygon(vector<vector<point2D> > obstacles, point2D genericPoint) {

	vector<point2D> visiblePolygon;
	vector<vector<lineSegment2D> > allPolygonEdges;

	for (int i = 0; i < obstacles.size(); i++) {
		vector<lineSegment2D> individualPolygonEdges = makePolygonEdges(obstacles[i]);
		allPolygonEdges.push_back(individualPolygonEdges); 
	} //Returns vector of all obstacle's edges.

	vector<eventPoint> intersectionEvents, tempContainer;

	lineSegment2D guardVertexEdge;
	eventPoint newEvent;

	int prevPointIndex, nextPointIndex;

	for (int i = 0; i < obstacles.size(); i++) { //For all obstacles.

		vector<point2D> individualPolygonPoints = obstacles[i];

		for (int j = 0; j < individualPolygonPoints.size(); j++) { //For all points inside an individual polygon

			if (j == 0) {
				prevPointIndex = individualPolygonPoints.size() - 1;
				nextPointIndex = 1;
			}

			else if (j == individualPolygonPoints.size() - 1) {
				prevPointIndex = individualPolygonPoints.size() - 2;
				nextPointIndex = 0;
			}

			else {
				prevPointIndex = j - 1;
				nextPointIndex = j + 1;
			}

			if (pointComparer(genericPoint, individualPolygonPoints[j])) {
				continue;
			}

			if (individualPolygonPoints[j].x == 0 && individualPolygonPoints[j].y == 0) {
				continue;
			}

			guardVertexEdge = extendLine(individualPolygonPoints[j], genericPoint);

			for (int k = 0; k < allPolygonEdges.size(); k++) { //For all polygon's edges vector of vector of lines

				vector<lineSegment2D> individualPolygonEdges = allPolygonEdges[k];

				for (int f = 0; f < individualPolygonEdges.size(); f++) { //For all edges in an individual polgyon

					if (linesIntersect(guardVertexEdge, individualPolygonEdges[f])) {
			 			newEvent.point = intersectionPoint(guardVertexEdge.p1, guardVertexEdge.p2, individualPolygonEdges[f].p1, individualPolygonEdges[f].p2);

				 		if (eventExists(intersectionEvents, newEvent.point)) {
							continue;
						}

					initializeNewEvent(&newEvent, f, individualPolygonEdges[f], guardVertexEdge, genericPoint);
			 		intersectionEvents.push_back(newEvent);

			 		}

				}
			}

			sort(intersectionEvents.begin(), intersectionEvents.end(), eventSorterDistance);

			for (int f = 0; f < intersectionEvents.size(); f++) {

				initializePoint(&intersectionEvents[f], genericPoint, individualPolygonPoints[prevPointIndex]);

				if (intersectionEvents[f].intersectionType == 1) {
					break;
				}
				
				else if (intersectionEvents[f].intersectionType == 0 && !pointComparer(intersectionEvents[f].point, genericPoint)) {
					visiblePolygon.push_back(intersectionEvents[f].point);
					break;
				}

			}

			intersectionEvents.clear();

		}

	}

	return visiblePolygon;
}

bool endVisible(point2D startPoint, point2D endPoint, vector<vector<point2D> > obstacles) {

	vector<point2D> visiblePolygon;
	vector<vector<lineSegment2D> > allPolygonEdges;

	for (int i = 0; i < obstacles.size(); i++) {
		vector<lineSegment2D> individualPolygonEdges = makePolygonEdges(obstacles[i]);
		allPolygonEdges.push_back(individualPolygonEdges); 
	} //Returns vector of all obstacle's edges.

	vector<eventPoint> intersectionEvents, tempContainer;

	lineSegment2D guardVertexEdge;
	eventPoint newEvent;

	guardVertexEdge = extendLine(startPoint, endPoint);

	for (int k = 0; k < allPolygonEdges.size(); k++) { //For all polygon's edges vector of vector of lines

		vector<lineSegment2D> individualPolygonEdges = allPolygonEdges[k];

		for (int f = 0; f < individualPolygonEdges.size(); f++) { //For all edges in an individual polgyon

			if (properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, individualPolygonEdges[f].p1, individualPolygonEdges[f].p2)) {
				return false;
			}
		}
	}

	return true;

}

void visibilityGraph (vector<vector<point2D> > *obstacles, point2D *startPoint, point2D *endPoint) {

	vector<point2D> visiblePoints, validPoints;
	vector<vector<lineSegment2D> > allPolygonEdges;
	vector<vector<point2D> > obstaclesVector = (*obstacles);

	vector<point2D> endP;
	endP.push_back((*endPoint));

	for (int i = 0; i < (*obstacles).size(); i++) {

		for (int j = 0; j < (*obstacles)[i].size(); j++) {

			visiblePoints = visiblePolygon((*obstacles), (*obstacles)[i][j]);

			for (int k = 0; k < visiblePoints.size(); k++) {
				(*obstacles)[i][j].visiblePoints.push_back(visiblePoints[k]);
			}

			if (endVisible((*obstacles)[i][j], (*endPoint), (*obstacles))) {
				(*obstacles)[i][j].visiblePoints.push_back((*endPoint));
			}

		}

	}

	visiblePoints = visiblePolygon((*obstacles), (*endPoint));

	for (int i = 0; i < visiblePoints.size(); i++) {
		(*endPoint).visiblePoints.push_back(visiblePoints[i]);
	}


	visiblePoints = visiblePolygon((*obstacles), (*startPoint));

	for (int i = 0; i < visiblePoints.size(); i++) {
		(*startPoint).visiblePoints.push_back(visiblePoints[i]);
	}

	if (endVisible((*startPoint), (*endPoint), (*obstacles))) {
		(*startPoint).visiblePoints.push_back((*endPoint));
	}

}

bool distanceSorter(const point2D& pointOne, const point2D& pointTwo) {
  return pointOne.distance < pointTwo.distance;

}

bool lineSegmentsIdentical(lineSegment2D l1, lineSegment2D l2) { //tested works

	if((l1.p1.x == l2.p1.x && l1.p1.y == l2.p1.y && l1.p2.x == l2.p2.x && l1.p2.y == l2.p2.y) || 
	(l1.p1.x == l2.p2.x && l1.p1.y == l2.p2.y && l1.p2.x == l2.p1.x && l1.p2.y == l2.p1.y)) {
		return true;
	} else {
		return false;
	}
}

bool operator<(const point2D& lhs, const point2D& rhs) {
	return lhs.distance > rhs.distance;  //<?
}

int previousIndex(point2D visiblePoint, graph myGraph) {

	for (int i = 0; i < myGraph.vertices.size(); i++) {

		if(pointComparer(visiblePoint, myGraph.vertices[i])) {
			return myGraph.vertices[i].global_index;
		}

	}

	return 0;

}

bool notVisited(point2D newPoint, vector<point2D> visitedPoints) {

	for (int i = 0; i < visitedPoints.size(); i++) {

		if (pointComparer(newPoint, visitedPoints[i])) {
			return false;
		}

	}

	return true;

}

vector<point2D*> dijkstra2(vector<point2D> obstaclePoints, point2D startPoint, point2D endPoint) {
	graph myGraph;
	myGraph.vertices = obstaclePoints;

	point2D temporaryPointer;

	vector<point2D> visitedPoints;

	startPoint.distance = 0;

	vector<point2D*> prevVector (myGraph.vertices.size(), NULL);

	vector<point2D> tempVector;

	priority_queue<point2D, vector<point2D> > pq;

	for(int i=0; i<myGraph.vertices.size(); ++i) {

		cout << "Index Point: " << "Point: X: " << myGraph.vertices[i].x << ", Y: " << myGraph.vertices[i].y << endl;
		cout << "I: " << i << endl;

		myGraph.vertices[i].global_index = i;
	}

	//cout << "Vertex [x: " << myGraph.vertices[0].x << " , y: " << myGraph.vertices[0].y << "]" << endl;

	//cout << "prevVector size: " << prevVector.size() << endl;

	for(int i=0; i<myGraph.vertices.size(); ++i) {
		for(int j=0; j<myGraph.vertices[i].visiblePoints.size(); ++j) {
			myGraph.vertices[i].visiblePoints[j].distance = 1000.0;
		}
	}

	pq.push(myGraph.vertices[0]); //i.e. start point

	// while (!pq.empty()) {
	// 	point2D tempPoint = pq.top();
	// 	pq.pop();
	// 	tempVector.push_back(tempPoint);
	// }

	// cout << "Points In Queue: " << endl;

	// for (int i = 0; i < tempVector.size(); i++) {

	// 	cout << "Point: X: " << tempVector[i].x << ", Y: " << tempVector[i].y << endl;

	// }

	while (!pq.empty()) {

		point2D minDistPoint = pq.top();
		pq.pop();

		visitedPoints.push_back(minDistPoint);


		// for (int i = 0; i < minDistPoint.visiblePoints.size(); i++) {

		// 	cout << "Visible Point: X: " << minDistPoint.visiblePoints[i].x << ", Y: " << minDistPoint.visiblePoints[i].y << endl;

		// }

		if (pointComparer(minDistPoint, startPoint)) {
				cout << "FIRST Point!" << endl;
			}

		else {

			cout << "NEW POINT" << endl;
		}

		printf("MinDistPoint ( %i , %i )\n",minDistPoint.x,minDistPoint.y);

		cout << "MinDistPoint Visible POints Size: " << minDistPoint.visiblePoints.size() << endl;

		for (int j = 0; j < minDistPoint.visiblePoints.size(); j++) {

			if (!notVisited(minDistPoint.visiblePoints[j], visitedPoints)) {
				continue;
			}

			for (int i = 0; i < myGraph.vertices.size(); i++) {

				if (pointComparer(myGraph.vertices[i], minDistPoint.visiblePoints[j])) {
					temporaryPointer = myGraph.vertices[i];
				}

			}

			//cout<<"VISITING NBORS\n";
			double alt = minDistPoint.distance + segmentLength(temporaryPointer, minDistPoint);
			//cout<<"GOT ALT\n";

			printf("MinDistPoint Neighbour( %i , %i )\n",minDistPoint.visiblePoints[j].x,minDistPoint.visiblePoints[j].y);


			// cout << "Segment Length: " << segmentLength(minDistPoint, minDistPoint.visiblePoints[j]) << endl;
			// cout << "Min Point Distance: " << minDistPoint.distance << endl;

			//double alt = segmentLength(minDistPoint, minDistPoint.visiblePoints[i]) + minDistPoint.distance;

			cout << "Alternate Distance: "<<alt<<endl;
			cout << "Distance to Point: "<< temporaryPointer.distance<<endl;



			if (alt < temporaryPointer.distance) {
				//cout<<"CHECK PASSED\n";
				minDistPoint.visiblePoints[j].distance = alt; //Need to change globally.

				cout << "New Distance: " << minDistPoint.visiblePoints[j].distance << endl;

				//cout << "Current Neighbour: "

				for (int i = 0; i < myGraph.vertices.size(); i++) { //Update distance in myGraph.vertices.
					if (pointComparer(myGraph.vertices[i], minDistPoint.visiblePoints[j])) {

						cout << "Vertex To Be Updated: " << "Visible Point: X: " << minDistPoint.visiblePoints[j].x << ", Y: " << minDistPoint.visiblePoints[j].y << endl;

						myGraph.vertices[i].distance = alt;	

						cout << "New Distance 2: " << myGraph.vertices[i].distance << endl;

					}
				}

				int prevIndex = previousIndex(minDistPoint.visiblePoints[j], myGraph);

				//cout << "Memory Space: " << prevVector[prevIndex] << endl;

				if(prevVector[prevIndex] == NULL) {

					cout << endl << endl << endl;

					cout << "PREv Index: " << prevIndex << endl;
					cout << "POINT TO : << Visible Point: X: " << minDistPoint.visiblePoints[j].x << ", Y: " << minDistPoint.visiblePoints[j].y << endl;
					cout << "POINT FROM : << Visible Point: X: " << minDistPoint.x << ", Y: " << minDistPoint.y << endl;

					cout << endl << endl << endl;

					for (int i = 0; i < myGraph.vertices.size(); i++) {

						if (pointComparer(myGraph.vertices[i], minDistPoint)) {
							prevVector[prevIndex] = &myGraph.vertices[i];

							cout << "PREV POINT. X: " << (*prevVector[prevIndex]).x << ", Y: " << (*prevVector[prevIndex]).y << endl;


						}

					}

					//prevVector[prevIndex] = &minDistPoint;

				}

				// if (pointComparer(minDistPoint.visiblePoints[j], endPoint)) {



				// 	return prevVector;
				// }

				else {
					cout << "Not found" << endl;
				}

				if (!pq.empty()) {
					while (!pq.empty()) {
						point2D tempPoint = pq.top();
						pq.pop();
						tempVector.push_back(tempPoint);
					}

				}

				else {

					cout << "empty!" << endl;

					cout << "Point To Be Added: " << "Point: X: " << minDistPoint.visiblePoints[j].x << ", Y: " << minDistPoint.visiblePoints[j].y << endl;

					for (int i = 0; i < myGraph.vertices.size(); i++) {

						if (pointComparer(myGraph.vertices[i],minDistPoint.visiblePoints[j])) {

							if (notVisited(myGraph.vertices[i], visitedPoints)) {
								tempVector.push_back(myGraph.vertices[i]);
								break;

							}

						}


					}

				}

				cout << "Points In Queue: " << endl;

				cout << tempVector.size() << endl;

				for (int i = 0; i < tempVector.size(); i++) {

					cout << "Point: X: " << tempVector[i].x << ", Y: " << tempVector[i].y << endl;

				}

				bool inQueue = false;

				for (int i = 0; i < tempVector.size(); i++) {
					if (pointComparer(tempVector[i], minDistPoint.visiblePoints[j])) {
						inQueue = true;

						break;

					}

				}

				if (inQueue == false) {

					for (int i = 0; i < myGraph.vertices.size(); i++) {

						if (pointComparer(myGraph.vertices[i],minDistPoint.visiblePoints[j])) {

							if (notVisited(myGraph.vertices[i], visitedPoints)) {

								cout << "ADDING NEW POINT" << endl;

									tempVector.push_back(myGraph.vertices[i]);
									break;

								}

						}
					}
				}

				for (int k = 0; k < tempVector.size(); k++) {

					cout << "Adding POint to Queue: " << "Point: X: " << tempVector[k].x << ", Y: " << tempVector[k].y << endl;

					pq.push(tempVector[k]);

				}
				tempVector.clear();

			}

	}

	}

	cout << "IN GEOM" << endl;

	for (int i = 0; i < prevVector.size(); i++) {

        cout << "I: " << i << endl;

        if (prevVector[i] != NULL) {

        	if ((*prevVector[i]).y > 5000) {

        		cout << "Start Point. X: " << startPoint.x << ", Y: " << startPoint.y << endl;

        		(*prevVector[i]).x = startPoint.x;
        		(*prevVector[i]).y = startPoint.y;
        	}

          cout << "Parent. X: " << (*prevVector[i]).x << ", Y: " << (*prevVector[i]).y << endl;

        }



      }



	return prevVector;


}