/* 
   Laura Toma
*/

#ifndef __geom_h
#define __geom_h

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <vector>
#include <map>
#include <utility>

using namespace std;

typedef struct _point2d {
  int x,y;
  double guardAngle;
  double distanceFromGuard; 
  bool left; // true/1 = left, false/0 = right. 
  int intersectionType; // 1 = proper, 0 = improper
  int intersectionEdge;
  int distanceFromPreviousPoint;
} point2D;

typedef struct _eventPoint {
	point2D point;
	int intersectionType; // 1 = proper, 0 = improper
	double distanceFromGuard;
  int intersectionEdge;
  int distanceFromPreviousPoint;
} eventPoint;

typedef struct _lineSegment2D {
    point2D p1, p2;
} lineSegment2D;


typedef struct _rect2D  {
    point2D origin;
    float width, height;
} rect2D;

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

bool isSimple(vector<point2D> polygonPoints);
vector<point2D> polygonShift(vector<point2D> originalPolygonPoints, point2D guardPoint);
double computeXIntersection (int firstPointX, int firstPointY, int secondPointX, int secondPointY);
bool isInPolygon(vector<point2D> polygonPoints, point2D guardPoint);
bool collinear(point2D a, point2D b, point2D c);
int Area2(point2D a, point2D b, point2D c);
bool Xor(bool x, bool y);
bool left(point2D a, point2D b, point2D c);
bool properIntersect(point2D a, point2D b, point2D c, point2D d);
bool between(point2D a, point2D b, point2D c);
bool improperIntersect(point2D a, point2D b, point2D c, point2D d);
double segmentLength(point2D point1, point2D point2);
lineSegment2D extendLine(point2D polygonPoint, point2D guardPoint);
lineSegment2D extendLineBothDirections(point2D polygonPoint, point2D guardPoint);
bool linesIntersect(lineSegment2D guardVertexEdge, lineSegment2D polygonEdge);
bool grazePoint(point2D polygonPointV1, point2D polygonPointV2, lineSegment2D guardVertexEdge);
point2D intersectionPoint(point2D a, point2D b, point2D c, point2D d);
bool eventSorterDistance(const eventPoint& pointOne, const eventPoint& pointTwo);
double pointDistance(point2D event, point2D guardPoint);
bool pointExists(vector<point2D> visiblePolygon, point2D newPoint);
bool eventExists(vector<eventPoint> events, point2D newPoint);
vector<lineSegment2D> makePolygonEdges(vector<point2D> polygonPoints);
bool pointComparer(const point2D& pointOne, const point2D& pointTwo);
bool distanceSort(const point2D& pointOne, const point2D& pointTwo);
void remove_coincident_points(vector<point2D> *points);
vector<point2D> reformPolygon(vector<point2D> visiblePolygonPoints, vector<lineSegment2D> polygonEdges);
void initializeNewEvent(eventPoint *newEvent, int intersectionEdge, lineSegment2D polygonEdge, lineSegment2D guardVertexEdge, point2D guardPoint);
void initializePoint(eventPoint *newEvent, point2D guardPoint, point2D prevPolygonPoint);
vector<point2D> visiblePolygon(vector<point2D> polygonPoints, point2D guardPoint);
vector<vector<point2D> > triangulatedPolygon(vector<point2D> polygonPoints, point2D guardPoint);

#endif
