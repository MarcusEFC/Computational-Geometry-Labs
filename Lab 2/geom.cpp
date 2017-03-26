#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <vector>
#include <map>
#include <utility>

using namespace std; 

//This multimap is the active structure containing the active horizontal segments
multimap<int, segment2D> activeStructure; 

/* isHorizontal
 *
 * Inputs:
 *     A line segment whose orientation need to be determined.
 *
 * Return Values:
 *     True if the segment is horizontal, false otherwise.
 *
 * Functional Description:
 *     For a line segment to be horizontal, the start point and end point y-coordinates need to be the same.
 *     This function checks to see if the two y-coordinates are equal and returns accordingly.
 */

int isHorizontal(segment2D checkedSegment) {

	if (checkedSegment.start.y == checkedSegment.end.y) {
		return 1;

	}

	else {
		return 0;

	}

}

/* getEvents
 * 
 * Inputs:
 *     - A vector of segments containing all line segments
 *     - A pointer to a vector of sweepEvents which this function will be adding events to.
 *
 * Return Values:
 *     N/A.
 *
 * Functional Description:
 *     This function iterates through all segments in the segments vector, and for each vector, determines its
 *     orientation using the isHorinzontal function, and continues accordingly. If the segment is horizontal,
 *     this function creates two new events, a horizontal segments start and end event, which includes its type
 *     (whether it is a start or end event), its x and y coordinates, and its associated segment, and inserts 
 *	   both into the sweep event vector. If the segment is vertical, it creates a single sweep event, storing
 *	   its type, x coordinate, and associated segment.
 */

void getEvents(vector<segment2D> segments, vector<sweepEvent> *events) {

	for (int i = 0; i < segments.size(); i++) {

		if (isHorizontal(segments[i])) {

			sweepEvent horizontalEventStart;
			sweepEvent horizontalEventEnd;

			horizontalEventStart.type = 0;
			horizontalEventStart.xCor = segments[i].start.x;
			horizontalEventStart.yCor = segments[i].start.y;
			horizontalEventStart.eventSegment = segments[i];

			horizontalEventEnd.type = 1;
			horizontalEventEnd.xCor = segments[i].end.x;
			horizontalEventEnd.yCor = segments[i].end.y;
			horizontalEventEnd.eventSegment = segments[i];

			(*events).push_back(horizontalEventStart);
			(*events).push_back(horizontalEventEnd);


		}

		else {
			sweepEvent vertical;

			vertical.type = 2;
			vertical.xCor = segments[i].start.x; // For vertical points, start.x == end.x
			vertical.yCor = -1; //y-coordinate does not need to be stored
			vertical.eventSegment = segments[i];

			(*events).push_back(vertical);

		}

	}

}

/* eventSorter
 * 
 * Inputs:
 *     Two points in the sweepEvent vector.
 *
 * Return Values:
 *     This function returns an int value representing which event has a smaller x coordinate.
 *
 * Functional Description:
 *     This function gives a description to the predefined sort function,
 *     as to how to sort the vector of points. We are sorting the vector
 *     based on the events x-coordinates.
 */

int eventSorter(const sweepEvent& eventOne, const sweepEvent& eventTwo) {

  return eventOne.xCor < eventTwo.xCor;

}

/* rangeSearch
 * 
 * Inputs:
 *     - The vertical line whose intersections needs to be determined.
 *	   - A pointer to the vector of points containing all intersection points that this function will populate.
 *
 * Return Values:
 *     N/A.
 *
 * Functional Description:
 *     This function first determines the lower and upper bound of the vertical segment in question. This function
 *     then iterates through the active structure multimap (the multimap containing all horizontal segments who are
 *     are currently being traversed) between the calculated upper and lower bound, and creates a new itersection for 
 *	   every horizontal segment currently in the active structure. This new intersection is inserted into the vector
 *     of intersection points.
 */

void rangeSearch(sweepEvent verticalLine, vector<point2D> *intersectionPoints) {

	multimap<int, segment2D>::iterator it, itLower, itUpper;
	int lowerY;
	int upperY;

	if (verticalLine.eventSegment.start.y > verticalLine.eventSegment.end.y) {
		upperY = verticalLine.eventSegment.start.y;
		lowerY = verticalLine.eventSegment.end.y;

	}

	else {
		upperY = verticalLine.eventSegment.end.y;
		lowerY = verticalLine.eventSegment.start.y;

	}

	itLower = activeStructure.lower_bound(lowerY);
	itUpper = activeStructure.upper_bound(upperY);

	for (it = itLower; it != itUpper; ++it) {

		point2D newIntersection;

		newIntersection.x = verticalLine.xCor;
		newIntersection.y = (*it).first;

		(*intersectionPoints).push_back(newIntersection);

	}

}

/* operator==
 * 
 * Inputs:
 *     Two segments which need ot be compared.
 *
 * Return Values:
 *     This function returns whether or not the two segments are the same.
 *
 * Functional Description:
 *     This function compares all the elements of the two segment structs, and if all are the same, then the two segments are the same, and the function returns true.
 *	   This function returns false otherwise.
 */

bool operator==(const segment2D& segmentOne, const segment2D& segmentTwo) {
    if (segmentOne.start.x == segmentTwo.start.x && segmentOne.end.x == segmentTwo.end.x && segmentOne.end.x == segmentTwo.end.x && segmentOne.end.y == segmentTwo.end.y) {
    	return true;
    }

    else {
    	return false;
    }

}

/* line_segment_scan
 * 
 * Inputs:
 *     - The vector containing all segments
 *	   - The vector containing all events which will be iterated over.
 *
 * Return Values:
 *     This function returns a vector of points containing all the intersection points.
 *
 * Functional Description:
 *     This function first calls the getEvents function, which populates the events vector with all the events which
 * 	   are to be iterated over. It then sorts the events based on their x-coordinate. Using a for loop, this function
 *	   iterates over all the events in the events vector, and based on the type of event (start, end or vertical segment),
 *     responds accordingly. If the event is the start of a horizontal segment, it inserts the segment into the active
 *	   structure with its y-coordinate as the key. If the events is the end of a horizontal segment, it removes the
 *     corresponding segment from the active structure, by iterating through all events in the map with the end event y
 *     coordinate, and using the operator== function to compare segments until the corresponding segment is found. When the
 *	   corrisponding event is found, that event is removed from the multimap. If the event is a vertical segment, the 
 *	   function calls rangeSearch, which populates the intersectionPoints vector.
 */

vector<point2D> line_segment_scan(vector<segment2D>  segments, vector<sweepEvent> events) {

	multimap<int, segment2D>::iterator it;

	vector<point2D> intersectionPoints;

	getEvents(segments, &events);

	sort(events.begin(), events.end(), eventSorter);

	for (int i = 0; i < events.size(); i++) {

		if (events[i].type == 0) { //Start event
			activeStructure.insert(pair<int, segment2D>(events[i].yCor, events[i].eventSegment));

		}

		else if (events[i].type == 1) { //End event

			int eventY = events[i].yCor;
			segment2D endSegment = events[i].eventSegment;

			multimap<int, segment2D>::iterator mapIt;
			pair <multimap<int, segment2D>::iterator, multimap<int, segment2D>::iterator> it;
			it = activeStructure.equal_range(eventY);

			for (mapIt = it.first; mapIt != it.second; ++mapIt) {

				segment2D mapSegment = mapIt->second;

				if (mapSegment == endSegment) {
					activeStructure.erase(mapIt);
					break;

				}

			}
			

		}

		else { //Vertical
			rangeSearch(events[i], &intersectionPoints);

		}


	}

	return intersectionPoints;

}
