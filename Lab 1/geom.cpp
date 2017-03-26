
#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stack>
#include <math.h>

using namespace std;

/* Inputs:
 *   Points a, b  and c which make up the triangle.
 *
 * Return Values:
 *     This function returns 2 times the area of the triangle with vertices a, b and c.
 *
 * Functional Description:
 *    Returns the signed area of triangle abc. The area is positive if c
 *    is to the left of ab, and negative if c is to the right of ab
 */

int signed_area2D(point2D a, point2D b, point2D c) { 
  return ((int)((a.x * (b.y - c.y)) - (a.y * (b.x - c.x)) + ((b.x * c.y) - (c.x * b.y))));
    
}

/* Inputs:
 *     Points p, q and r, for which you want to determine their collinearity.
 *
 * Return Values:
 *     This function returns an int specifying whether or not the points are collinear.
 *
 * Functional Description:
 *     Uses the property that, if the area of the triangle made up of points p, q, and r
 *     is zero, the points are collinear. This function returns accordingly.
 */

int collinear(point2D p, point2D q, point2D r) {
  return (signed_area2D(p, q, r) == 0);
    
}

/* Inputs:
 *     Points a, b, and c.
 *
 * Return Values:
 *     This function returns an int specifying whether or not point c is to the left of points a and b.
 *
 * Functional Description:
 *     Return 1 if c is strictly left of ab; 0 otherwise.
 */

int left (point2D a, point2D b, point2D c) {
  return ((((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y))) > 0);
    
}

/* Inputs:
 *     The vector containing all points.
 *
 * Return Values:
 *     This function returns the the point with the lowest y-coordinate to be the starting point.
 *
 * Functional Description:
 *     Iteratively checks points to find point with lowest y-coordinate. 
 *     If y-coordinates are equal, it returns the point with the lowest x-coordinate.
 */

point2D findFirstPoint(vector<point2D> pointVector) {

  point2D firstPoint = pointVector[0];

  for (int i = 1; i < pointVector.size(); ++i) {
    if (pointVector[i].y < firstPoint.y) {
      firstPoint = pointVector[i];

    }

    else if (pointVector[i].y == firstPoint.y) {

      if (pointVector[i].x < firstPoint.x) {
        firstPoint = pointVector[i];
      }

    }

  }

  return firstPoint;

}

/* Inputs:
 *     The reference point, and the point for which the ccw angle will be calculated.
 *
 * Return Values:
 *     This function returns the angle between the reference point and the new point.
 *
 * Functional Description:
 *     By computing the arc tangent between the two points, we return the angle between them.
 */

double computedAngle(point2D referencePoint, point2D newPoint) {
  return atan2((newPoint.y - referencePoint.y), (newPoint.x - referencePoint.x));

}

/* Inputs:
 *     The reference point, and the point for which the ccw angle will be calculated.
 *
 * Return Values:
 *     This function returns the computed distance between two points.
 *
 * Functional Description:
 *     By utilizing the distance formula equation of a line, we can calculate the distance between two points.
 */

double computedDistance(point2D referencePoint, point2D newPoint) {
  return sqrt((pow((newPoint.x - referencePoint.x), 2)) + (pow((newPoint.y - referencePoint.y), 2)));

}

/* Inputs:
 *     A pointer to the vector containing all points.
 *
 * Return Values:
 *     N/A. As we are using a pointer to the vector, we are modifying the vector directly.
 *
 * Functional Description:
 *     Calculates the angle and distance of the point from the reference point, and stores
 *     this information in the point struct.
 */

void computeProperties(vector<point2D> *pointVector) {

  point2D firstPoint = findFirstPoint((*pointVector));

  for (int i = 0; i < (*pointVector).size(); ++i) {

    if ((*pointVector)[i].x == firstPoint.x && (*pointVector)[i].y == firstPoint.y) {
      (*pointVector)[i].angle = -1;

    }

    else {
      double computeAngle = computedAngle(firstPoint, (*pointVector)[i]);

      if (computeAngle < 0) {
        computeAngle = computeAngle + 2*M_PI;

      }

      (*pointVector)[i].angle = computeAngle;

    }

    (*pointVector)[i].distance = computedDistance(firstPoint,(*pointVector)[i]);

  }

}

/* Inputs:
 *     Two points in the vector.
 *
 * Return Values:
 *     This function returns a boolean value representing whether or not 
 *     one point has a lower distance/angle than another point.
 *
 * Functional Description:
 *     This function gives a description to the predefined sort function,
 *     as to how to sort the vector of points. We are sorting the vector
 *     first based on the angle wrt to the reference point, and then the
 *     distance wrt to the reference point.
 */

bool pointSorter(const point2D& pointOne, const point2D& pointTwo) {
    
  if (pointOne.angle == pointTwo.angle) {
    return pointOne.distance < pointTwo.distance;

  }

  return pointOne.angle < pointTwo.angle;

}

/* Inputs:
 *     Pointer to the stack, and a pointer to the vector that encapsulates all of the points.
 *
 * Return Values:
 *     N/A.
 *
 * Functional Description:
 *     Removes a point struct from the vector, and pushes it onto the stack.
 */

void pushToStack(stack<point2D> *stackPointer, vector<point2D> *vectorPointer) {
  point2D auxPoint;
  
  auxPoint = (*vectorPointer).front();
  (*vectorPointer).erase((*vectorPointer).begin());
  (*stackPointer).push(auxPoint);
    
}

/* Inputs:
 *     The stack containing all of the points on the hull.
 *
 * Return Values:
 *     Returns the point below the point on the top of the stack.
 *
 * Functional Description:
 *     Pops the first point off the stack and returns the next point.
 */

point2D firstPoint(stack<point2D> hullStack) {
  hullStack.pop();
  point2D pointOne = hullStack.top();

  return pointOne;

}

/* Inputs:
 *     Vector p which contains all of the points generated.
 *
 * Return Values:
 *     Returns all of the points in the hull.
 *
 * Functional Description:
 *     First sorts the points and initializes the reference point to be the point with the lowest y-coordinate.
 *     Then initializes the stack and pushes points onto the stack in a counter-clockwise fashion, checking if
 *     the next point is to the left of the current two top items on the stack; if so, it pushes the item on the
 *     stack (this point may be apart of the hull), if not, i.e. the point is to the right, pop the top point in 
 *     the stack until the following point in the vector is to the left of the top two points. Repeat until all
 *     points in the vector have been checked. This function also takes into account that points can be collinear,
 *     where if three points are collinear, this function pops off the current top item, and pushes the next item
 *     onto the stack.
 */

vector<point2D> graham_scan(vector<point2D> p) {
  vector<point2D> result;
  stack<point2D> hullStack;
  point2D pointOne;

  computeProperties(&p);

  sort(p.begin(), p.end(), pointSorter);

  if (p.size() < 3) {
    return p;
  
  }

  else {
    pushToStack(&hullStack, &p);
    pushToStack(&hullStack, &p);

    while (p.size() != 0) {
        pointOne = firstPoint(hullStack);

        if (collinear(pointOne, hullStack.top(), p.front())) {
          hullStack.pop();
          pushToStack(&hullStack, &p);

        }
        
        else if (left(pointOne, hullStack.top(), p.front())) {
          pushToStack(&hullStack, &p);
            
        }
        
        else {

          while (!left(pointOne, hullStack.top(), p.front())) {
            hullStack.pop();
            pointOne = firstPoint(hullStack);

          }
          
          pushToStack(&hullStack, &p); 
            
        }  
        
    }

    while (hullStack.size() != 0) {
        point2D auxPoint = hullStack.top();
        hullStack.pop();
        result.push_back(auxPoint);

    }

  }

  return result; 
    
}
