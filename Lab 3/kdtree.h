
/* 
   LT based on code from John Visentin 
*/ 

#ifndef __kdtree_h
#define __kdtree_h


#include "geom.h"
#include <stddef.h> //for size_t
#include <vector>

using namespace std;

typedef struct _treeNode treeNode;

struct _treeNode {
  point2D p; /* If this is a leaf node, p represents the point
		stored in this leaf.  If this is not a leaf node, p
		represents the horizontal or vertical line stored in
		this node. For a vertical line, p.y is ignored. For
		a horizontal line, p.x is ignored
	     */
  char type; /* this can be 'h' (horizontal) or 'v' (vertical), or
		'l' (leaf) depending whether the node splits with
		a horizontal line or vertical line.  Technically
		this should be an enum.
	     */
  treeNode  *left, *right; /* left/below and right/above children. */
  
  /* we might need a parent pointer but we'll cross that bridge when
     we get to it */
};


typedef struct _kdtree {
  treeNode *root;  //root of the tree
  size_t count; //number of nodes in the tree
  size_t height; //height of the tree
} kdtree;

/* returns the root node of the given tree */
treeNode* kdtree_getRoot(kdtree *tree);

/* returns the point stored by the node */
point2D treeNode_getPoint(treeNode *node);

/* create a new empty tree */
kdtree* kdtree_init();

void treeNode_print(treeNode* node);

void treeNode_free(treeNode* node);

/* print out information about the tree including height, number of
 nodes, and each node in an in-order traversal */
void kdtree_print(kdtree *tree);

/* free all memory allocated for the tree, including the tree
   itself */
void kdtree_free(kdtree *tree);

int pointSorter(const point2D& pointOne, const point2D& pointTwo);

bool pointComparer(const point2D& pointOne, const point2D& pointTwo);

int pointSorterX(const point2D& pointOne, const point2D& pointTwo);

int pointSorterY(const point2D& pointOne, const point2D& pointTwo);

void remove_coincident_points(vector<point2D> *points);

treeNode* initializeNode(int x, int y, int depth, treeNode  *left, treeNode *right);

int medianIndex (int vectorSize);

vector<point2D> reInitializeSubsetXCor(vector<point2D> *initialX);

vector<point2D> reInitializeSubsetYCor(vector<point2D> *initialY);

/* create a new tree representing the given array of points */
kdtree* kdtree_build(vector<point2D> points);

treeNode* kdtree_build_rec(vector<point2D> sortedX, vector<point2D> sortedY, int depth);

#endif
