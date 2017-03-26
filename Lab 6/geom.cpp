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

bool pointComparer(const point2D& pointOne, const point2D& pointTwo) {

  if (pointOne.x == pointTwo.x && pointOne.y == pointTwo.y) {
    return true;
  }

  return false;

}

void remove_coincident_points(vector<point2D> *points) {
  vector<point2D>::iterator unique_end = (unique((*points).begin(), (*points).end(), pointComparer));
  (*points).erase(unique_end, (*points).end());
  
}

bool samePolygon(point2D startPoint, point2D endPoint, vector<point2D> polygon) {

	bool startOn, endOn = false;

	for (int i = 0; i < polygon.size(); i++) {

		if (pointComparer(startPoint, polygon[i])) {
			startOn = true;
		}

	}

	for (int i = 0; i < polygon.size(); i++) {

		if (pointComparer(endPoint, polygon[i])) {
			endOn = true;
		}

	}

	if (startOn && endOn) {
		return true;
	}

	return false;

}

bool midPointInside (point2D startPoint, point2D endPoint, vector<point2D> polygon) {

	point2D midPoint;

	midPoint.x = (startPoint.x + endPoint.x) / 2;
	midPoint.y = (startPoint.y + endPoint.y) / 2;

	if (isInPolygon(polygon, midPoint)) {
		return true;
	}

	return false;

}

bool areNeighbours(point2D startPoint, point2D endPoint, vector<point2D> polygon) {

	int startIndex, endIndex;

	for (int i = 0; i < polygon.size(); i++) {

		if (pointComparer(startPoint, polygon[i])) {
			startIndex = i;
		}

		else if (pointComparer(endPoint, polygon[i])) {
			endIndex = i;
		}

	}

	if ((startIndex + 1 == endIndex) || (endIndex + 1 == startIndex)) {
		return true;
	}

	return false;

}

vector<point2D> visiblePolygon(vector<vector<point2D> > obstacles, point2D genericPoint) {

	vector<point2D> visiblePolygon;
	vector<vector<lineSegment2D> > allPolygonEdges;

	for (int i = 0; i < obstacles.size(); i++) {
		vector<lineSegment2D> individualPolygonEdges = makePolygonEdges(obstacles[i]);
		allPolygonEdges.push_back(individualPolygonEdges); 
	} //Returns vector of all obstacle's edges.

	lineSegment2D guardVertexEdge;

	for (int i = 0; i < obstacles.size(); i++) { //For all obstacles.

		vector<point2D> individualPolygonPoints = obstacles[i];

		for (int j = 0; j < individualPolygonPoints.size(); j++) { //For all points inside an individual polygon

			if (pointComparer(genericPoint, individualPolygonPoints[j])) {
				continue;
			}

			if (individualPolygonPoints[j].x == 0 && individualPolygonPoints[j].y == 0) {
				continue;
			}

			if (samePolygon(individualPolygonPoints[j], genericPoint, individualPolygonPoints) && (areNeighbours(individualPolygonPoints[j], genericPoint, individualPolygonPoints))) {
				visiblePolygon.push_back(individualPolygonPoints[j]);
				continue;
			}

			guardVertexEdge.p1 = individualPolygonPoints[j];
			guardVertexEdge.p2 = genericPoint;

			bool properIntersection = false;

			for (int k = 0; k < allPolygonEdges.size(); k++) { //For all polygon's edges vector of vector of lines

				vector<lineSegment2D> individualPolygonEdges = allPolygonEdges[k];

				for (int f = 0; f < individualPolygonEdges.size(); f++) { //For all edges in an individual polgyon

					if (properIntersect(individualPolygonEdges[f].p1, individualPolygonEdges[f].p2, guardVertexEdge.p1, guardVertexEdge.p2)) {
						properIntersection = true;
						break;
					}

				}

			}

			if (!properIntersection) {
				if (samePolygon(individualPolygonPoints[j], genericPoint, individualPolygonPoints)) {
					if (!areNeighbours(individualPolygonPoints[j], genericPoint, individualPolygonPoints)) {
						if (!midPointInside(individualPolygonPoints[j], genericPoint, individualPolygonPoints)) {
							visiblePolygon.push_back(individualPolygonPoints[j]);
						}
						else {
							continue;
						}
					}
					else {
						visiblePolygon.push_back(individualPolygonPoints[j]);
					}
				}
				else {
					visiblePolygon.push_back(individualPolygonPoints[j]);
				}
			}

		}

	}

	return visiblePolygon;
}

bool endVisible(point2D startPoint, point2D endPoint, vector<vector<point2D> > obstacles) {

	vector<vector<lineSegment2D> > allPolygonEdges;

	for (int i = 0; i < obstacles.size(); i++) {
		vector<lineSegment2D> individualPolygonEdges = makePolygonEdges(obstacles[i]);
		allPolygonEdges.push_back(individualPolygonEdges); 
	} //Returns vector of all obstacle's edges.

	for (int i = 0; i < allPolygonEdges.size(); i++) { //For all polygon's edges vector of vector of lines

		vector<lineSegment2D> individualPolygonEdges = allPolygonEdges[i];

		for (int j = 0; j < individualPolygonEdges.size(); j++) { //For all edges in an individual polgyon

			if (properIntersect(startPoint, endPoint, individualPolygonEdges[j].p1, individualPolygonEdges[j].p2)) {
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

bool operator<(const point2D& lhs, const point2D& rhs) {
	return lhs.distance > rhs.distance;  //<?
}

int pointIndex(point2D point, vector<point2D> obstaclePoints) {

	for (int i = 0; i < obstaclePoints.size(); i++) {
		if (pointComparer(point, obstaclePoints[i])) {
			return i;
		}
	}

	return 0;

}

bool alreadyInQueue(priority_queue<point2D, vector<point2D> > pq, point2D point) {

	point2D tempPoint;

	while (!pq.empty()) {
		tempPoint = pq.top();
		pq.pop();

		if (pointComparer(tempPoint, point)) {
			return true;
		}

	}

	return false;

}

vector<point2D> dijkstraPQ(vector<point2D> obstaclePoints) {
	
	priority_queue<point2D, vector<point2D> > pq;
	point2D tempPoint;
	vector<point2D> prevVector(obstaclePoints.size());

	pq.push(obstaclePoints[0]);

	while (!pq.empty()) {

		tempPoint = pq.top();
		pq.pop();

		for (int i = 0; i < tempPoint.visiblePoints.size(); i++) {

			int correspondingIndex = pointIndex(tempPoint.visiblePoints[i], obstaclePoints);
			double alt = tempPoint.distance + segmentLength(tempPoint, obstaclePoints[correspondingIndex]);

			if (alt < obstaclePoints[correspondingIndex].distance) {
				obstaclePoints[correspondingIndex].distance = alt;
				prevVector[correspondingIndex] = tempPoint;

				if (!alreadyInQueue(pq, obstaclePoints[correspondingIndex])) {
					pq.push(obstaclePoints[correspondingIndex]);

				}

			}

		}

	}

	return prevVector;

}

vector<point2D> shortestPathPoints(vector<point2D> obstaclePoints, point2D startPoint, point2D endPoint) {
	
	vector<point2D> bigVector, pathPoints, shortestPath;
	bigVector = obstaclePoints;
	bigVector.push_back(startPoint);
	bigVector.push_back(endPoint);

	for (int i = 0; i < bigVector.size(); i++) {

        if (pointComparer(bigVector[i], startPoint)) {
          bigVector[i].distance = 0;
        }

        else if (pointComparer(bigVector[i], endPoint)) {
          bigVector[i].distance = 2001;
        }

        else {
          bigVector[i].distance = 2000;
        }

    }

    sort(bigVector.begin(), bigVector.end(), distanceSorter);

    pathPoints = dijkstraPQ(bigVector);

	point2D tempPoint = endPoint;

	bool foundStart = false;
	int nextIndex;

	shortestPath.push_back(endPoint);

	while(!foundStart) {

	    if (pointComparer(tempPoint, endPoint)) { //endPoint
	      tempPoint = pathPoints[pathPoints.size() - 1];
	      shortestPath.push_back(tempPoint);

	    }

	    for (int i = 0; i < bigVector.size(); i++) {
	      if (pointComparer(tempPoint, bigVector[i])) {
	        nextIndex = i;

	      }
	    }

	    tempPoint = pathPoints[nextIndex];

	    shortestPath.push_back(tempPoint);

	    if (pointComparer(tempPoint, startPoint)) {
	      foundStart = true;
	    }
	}

	return shortestPath;

}