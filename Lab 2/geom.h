#ifndef __geom_h
#define __geom_h

#include <vector>
#include <map>

using namespace std; 

typedef struct _point2d {
  int x,y; 
} point2D;


typedef struct _segment2d {
  point2D start; 
  point2D end; 
} segment2D;

typedef struct _event {
	int type; // 0 == start, 1 == end, 2 == vertical
	int xCor; // x-coordinate
	int yCor; // y-coordinate
	segment2D eventSegment;

} sweepEvent;

int isHorizontal(segment2D checkedSegment);

void getEvents(vector<segment2D> segments, vector<sweepEvent> *events);

int eventSorter(const sweepEvent& eventOne, const sweepEvent& eventTwo);

vector<point2D> line_segment_scan(vector<segment2D>  segments, vector<sweepEvent> events);

void rangeSearch(sweepEvent verticalLine, vector<point2D> *intersectionPointsPointer);

#endif
