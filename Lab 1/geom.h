#ifndef __geom_h
#define __geom_h

#include <vector>
#include <stack>

using namespace std; 

typedef struct _point2d {
  int x,y;
  double angle; //The angle between the reference point and the point.
  double distance; //The distance between the reference point and the point.

 } point2D;

int signed_area2D(point2D a, point2D b, point2D c); 

int collinear(point2D p, point2D q, point2D r);

int left (point2D a, point2D b, point2D c);

point2D findFirstPoint(vector<point2D> pointVector); 

double computedAngle(point2D referencePoint, point2D newPoint);

double computedDistance(point2D referencePoint, point2D newPoint);

void computeProperties(vector<point2D> *pointVector);

bool pointSorter(const point2D& pointOne, const point2D& pointTwo);

void pushToStack(stack<point2D> *stackPointer, vector<point2D> *vectorPointer);

point2D firstPoint(stack<point2D> pointStack);

vector<point2D> graham_scan(vector<point2D>);

#endif
