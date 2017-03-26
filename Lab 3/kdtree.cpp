/* 

Marcus Christiansen, based on code by John Visentin

*/

#include "kdtree.h"
#include <cstdlib>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <new>
#include <set>

using namespace std;

/* returns the root node of the given tree */
treeNode* kdtree_getRoot(kdtree *tree) {
  assert(tree); 
  return tree->root;

}

/* returns the point stored by the node */
point2D treeNode_getPoint(treeNode *node) {
  assert(node); 
  return node->p;

}

/* create a new empty tree */
kdtree* kdtree_init() {
  
  kdtree* tree = new (nothrow) kdtree;
  assert(tree); 

  tree->root = NULL; 
  tree->count = tree->height = 0; 
  return tree; 

}

/* private function to recursively print a node and its subtree in
   order */
void treeNode_print(treeNode* node) {

  if (node == NULL) return; 

  //if we are here, node must be valid

  //recursively print left child 
  treeNode_print(node->left); 

  //print node 
  switch (node->type) {
  case 'h':
    printf("HORIZONTAL: (y=%d)\n", node->p.y);
    break; 
  case 'v':
    printf("VERTICAL: (x=%d)\n", node->p.x); 
    break; 
  case 'l': 
    printf("LEAF: (p=(%d,%d))\n", node->p.x, node->p.y); 
    break;
  default: 
    printf("Improper tree node type\n"); 
    exit(1); 
  };

  //recursively print right child
  treeNode_print(node->right); 

} 

/* print out information about the tree including height, number of
 nodes, and each node in an in-order traversal */
void kdtree_print(kdtree *tree) {

  if (tree) {
    printf("--- kdtree Info ---\n");
    printf("Height: %lu, Node Count: %lu\n", tree->height, tree->count);
    printf("Nodes in order:\n\n");
    treeNode_print(tree->root);
    printf("---------------------\n");
  }
}

//private function to recursively free the subtree rooted at node
void treeNode_free(treeNode* node) {

  treeNode* leftChild;
  treeNode* rightChild;

  if (node->left == NULL && node->right == NULL) {
    delete node;
  }

  else if (node->left == NULL) {
    rightChild = node->right;
    treeNode_free(rightChild);

  }

  else if (node->right == NULL) {
    leftChild = node->left;
    treeNode_free(leftChild);

  }

  else {
    leftChild = node->left;
    rightChild = node->right;

    treeNode_free(leftChild);
    treeNode_free(rightChild);

  }

}

/* free all memory allocated for the tree, including the tree
   itself */
void kdtree_free(kdtree *tree) {
  if (!tree) return;
  treeNode_free(tree->root); 
  delete tree; 
  
}

/* pointSorter
 * 
 * Inputs:
 *     Two points.
 *
 * Return Values:
 *     This function returns an int value representing which point has a smaller x coordinate. If they both have the same, this returns
 *     which point has the smalle y coordinate.
 *
 * Functional Description:
 *     This function gives a description to the predefined sort function,
 *     as to how to sort the vector of points. We are sorting the vector
 *     based on the events x-coordinates first, and then by y.
 */

int pointSorter(const point2D& pointOne, const point2D& pointTwo) {

  if (pointOne.x == pointTwo.x) {
    return pointOne.y < pointTwo.y;
  }

  else {
    return pointOne.x < pointTwo.x;
  }

}

/* pointSorterX
 * 
 * Inputs:
 *     Two points.
 *
 * Return Values:
 *     This function returns an int value representing which point has a smaller x coordinate. 
 *
 * Functional Description:
 *     This function gives a description to the predefined sort function,
 *     as to how to sort the vector of points. We are sorting the vector
 *     based on the events x-coordinates.
 */

int pointSorterX(const point2D& pointOne, const point2D& pointTwo) {
  return pointOne.x < pointTwo.x;

}

/* pointSorterY
 * 
 * Inputs:
 *     Two points.
 *
 * Return Values:
 *     This function returns an int value representing which point has a smaller y coordinate. 
 *
 * Functional Description:
 *     This function gives a description to the predefined sort function,
 *     as to how to sort the vector of points. We are sorting the vector
 *     based on the events y-coordinates.
 */

int pointSorterY(const point2D& pointOne, const point2D& pointTwo) {
  return pointOne.y < pointTwo.y;

}

/* pointComparer
 * 
 * Inputs:
 *     Two points which need ot be compared.
 *
 * Return Values:
 *     This function returns whether or not the two points are the same.
 *
 * Functional Description:
 *     This function compares all the x and y cooridnates of the two points, and if all are the same, then the two segments are the same, 
 *     and the function returns true. This function returns false otherwise.
 */

bool pointComparer(const point2D& pointOne, const point2D& pointTwo) {

  if (pointOne.x == pointTwo.x && pointOne.y == pointTwo.y) {
    return true;
  }

  else {
    return false;
  }

}

/* remove_coincident_points
 * 
 * Inputs:
 *     A pointer to a vector of points.
 *
 * Return Values:
 *     void.
 *
 * Functional Description:
 *     This function first sorts the vector of points being pointed to, and then removes all instances of identical points using the
 *     pre-defined 'unique' function, using our pointComparer function to determine if two points are coincident.
 */

void remove_coincident_points(vector<point2D> *points) {
  sort((*points).begin(), (*points).end(), pointSorter);
  vector<point2D>::iterator unique_end = (unique((*points).begin(), (*points).end(), pointComparer));
  (*points).erase(unique_end, (*points).end());
  
}

/* initializeNode
 * 
 * Inputs:
 *     - x coordinate of new node in tree.
 *     - y coordinate of new node in tree.
 *     - depth of current state of tree.
 *     - pointer to new node's left child.
 *     - pointer to new node's right child.
 *
 * Return Values:
 *     A treeNode pointer pointing to a new node in the tree.
 *
 * Functional Description:
 *     This function initializes a new node that is to be inserted into the tree, by creating a new treeNode pointer, initializing all of its properties using
 *     the function inputs, and returns a pointer to the newly created treeNode.
 */

treeNode* initializeNode(int x, int y, int depth, treeNode  *left, treeNode *right) { //CHANGE TO RETURN POINTER.

  treeNode* treeLeaf = new (nothrow) treeNode; //NEED TO MALLOC!!! //RETURN THE POINTER!!!

  assert(treeLeaf);

  point2D nodePoint;
  nodePoint.x = x;
  nodePoint.y = y;

  treeLeaf->p = nodePoint;

  if ((depth % 2) == 0) {
    treeLeaf->type = 'v';
  }

  else {
    treeLeaf->type = 'h';
  }

  treeLeaf->left = left;

  treeLeaf->right = right;

  return treeLeaf;

}

/* medianIndex
 * 
 * Inputs:
 *     A vector size (int)
 *
 * Return Values:
 *     This function returns the index of the median value in a vector depending on if the vector has an even or odd number of entries.
 *
 * Functional Description:
 *     This function determines if the vector size is even or odd, and returns the appropriate median index in the vector.
 */

int medianIndex (int vectorSize) {

  if (vectorSize % 2 == 0) {
    return (vectorSize / 2) - 1;
  }

  else {
    return vectorSize / 2;
  }

}

/* reInitializeSubsetXCor
 * 
 * Inputs:
 *     A pointer to a vector of points sorted by y coordinate.
 *
 * Return Values:
 *     A vector of points above and not including the median point with respect to the point's x coordinate.
 *
 * Functional Description:
 *     This function first calculates the median index of the vector being pointed, and based on this value, divides the array being pointed to into two subarrays,
 *     by creating a new vector of points, and storing all points below and including the median point, and removing all points below and including the median point
 *     in the vector being pointed to, resulting in a vector of all points above the median point. This way two sub vectors are created; one is returned, and the other
 *     is created by manipulating the vector being pointed to. 
 */

vector<point2D> reInitializeSubsetXCor(vector<point2D> *initialY) {

  vector<point2D> subYLo;

  sort((*initialY).begin(), (*initialY).end(), pointSorterX);

  int vectorSize = (*initialY).size();

  int medIndex = medianIndex(vectorSize);

  for (int i = 0; i <= medIndex; i ++) {
    subYLo.push_back((*initialY)[i]);
  }

  (*initialY).erase((*initialY).begin(), (*initialY).begin() + medIndex + 1);

  sort((*initialY).begin(), (*initialY).end(), pointSorterY);
  sort(subYLo.begin(), subYLo.end(), pointSorterY);

  return subYLo;

}

/* reInitializeSubsetXCor
 * 
 * Inputs:
 *     A pointer to a vector of points sorted by x coordinate.
 *
 * Return Values:
 *     A vector of points above and not including the median point with respect to the point's y coordinate.
 *
 * Functional Description:
 *     Same function as reInitializeSubsetXCor, however, this function divides the vector being pointed to based on the point's y value.
 */

vector<point2D> reInitializeSubsetYCor(vector<point2D> *initialX) {

  vector<point2D> subXLo;

  sort((*initialX).begin(), (*initialX).end(), pointSorterY);

  int vectorSize = (*initialX).size();

  int medIndex = medianIndex(vectorSize);

  for (int i = 0; i <= medIndex; i ++) {
      subXLo.push_back((*initialX)[i]);
  }

  (*initialX).erase((*initialX).begin(), (*initialX).begin() + medIndex + 1);

  sort((*initialX).begin(), (*initialX).end(), pointSorterX);
  sort(subXLo.begin(), subXLo.end(), pointSorterX);

  return subXLo;
  
}

/* kdtree_build_rec
 * 
 * Inputs:
 *     - A vector of all points sorted by their x coordinate
 *     - A vector of all points sorted by their y coordinate
 *     - The current depth of the tree
 *
 * Return Values:
 *     This function returns a pointer to the root of the kdtree.
 *
 * Functional Description:
 *     This function recursively constructs the kdtree based on the number of remaining points and the current depth of the tree. If the size of the vector
 *     is 1, then the node is a leaf node, and the function creates a new treeNode pointer, initializes the node being pointed to, and returns the new pointer.
 *     Otherwise, if the depth is even, this function creates a new vertical line using the median x-coordinate point (which is determined according to if there
 *     are an even or odd number of points), and recursively calls the build function twice, which will return the node's left and right child respectively. 
 *     Otherwise, if the depth is odd, this function creates a new horizontal line, and recursively calls the build function twice, which will again return 
 *     the node's left and right child respectively. This ultimately builds the kdtree for the series of points inputted.
 */

treeNode* kdtree_build_rec(vector<point2D> sortedX, vector<point2D> sortedY, int depth) {

  if (sortedX.size() == 1) { //Leaf Node

    treeNode* rootPointer  = new (nothrow) treeNode;
    rootPointer = initializeNode(sortedX[0].x, sortedY[0].y, depth, NULL, NULL);

    rootPointer->type = 'l';

    return rootPointer;

  }

  else if ((depth % 2) == 0) { //Vertical Line (median of x points)

    int treeDepth = depth;

    vector<point2D> subYLo, subYHi;

    if ((sortedX.size() % 2) == 0) { //Even number of points

      size_t const splitSize = medianIndex(sortedX.size()); // (vectorSize / 2) - 1;

      point2D internalPoint = sortedX[splitSize];

      vector<point2D> subXLo(sortedX.begin(), sortedX.begin() + splitSize + 1);
      vector<point2D> subXHi(sortedX.begin() + splitSize + 1, sortedX.end());
      subYLo = reInitializeSubsetXCor(&sortedY);
      subYHi = sortedY;

      treeNode* leftChild = kdtree_build_rec(subXLo, subYLo, treeDepth + 1);
      treeNode* rightChild = kdtree_build_rec(subXHi, subYHi, treeDepth + 1);

      treeNode *internalNode = initializeNode(internalPoint.x, internalPoint.y, depth, leftChild, rightChild);

      return internalNode;

    }

    else { //Odd number of points

      int treeDepth = depth;

      size_t const splitSize = medianIndex(sortedX.size()); // vectorSize / 2;

      point2D internalPoint = sortedX[splitSize];

      vector<point2D> subXLo(sortedX.begin(), sortedX.begin() + splitSize + 1);
      vector<point2D> subXHi(sortedX.begin() + splitSize + 1, sortedX.end());
      subYLo = reInitializeSubsetXCor(&sortedY);
      subYHi = sortedY;

      treeNode* leftChild = kdtree_build_rec(subXLo, subYLo, treeDepth + 1);
      treeNode* rightChild = kdtree_build_rec(subXHi, subYHi, treeDepth + 1);

      treeNode *internalNode = initializeNode(internalPoint.x, internalPoint.y, depth, leftChild, rightChild);

      return internalNode;

    }

  }

  else { //Horizontal Line (median of y points)

    int treeDepth = depth;

    vector<point2D> subXLo, subXHi;

    if ((sortedY.size() % 2) == 0) { //Even number of points

      size_t const splitSize = medianIndex(sortedY.size()); // (vectorSize / 2) - 1;

      point2D internalPoint = sortedY[splitSize];

      vector<point2D> subYLo(sortedY.begin(), sortedY.begin() + splitSize + 1); 
      vector<point2D> subYHi(sortedY.begin() + splitSize + 1, sortedY.end());
      subXLo = reInitializeSubsetYCor(&sortedX);
      subXHi = sortedX;

      treeNode* leftChild = kdtree_build_rec(subXLo, subYLo, treeDepth + 1);
      treeNode* rightChild = kdtree_build_rec(subXHi, subYHi, treeDepth + 1);

      treeNode *internalNode = initializeNode(internalPoint.x, internalPoint.y, depth, leftChild, rightChild);

      return internalNode;

    }

    else { //Odd number of points

      int treeDepth = depth;

      size_t const splitSize = medianIndex(sortedY.size()); // vectorSize / 2;

      point2D internalPoint = sortedY[splitSize];

      vector<point2D> subYLo(sortedY.begin(), sortedY.begin() + splitSize + 1); 
      vector<point2D> subYHi(sortedY.begin() + splitSize + 1, sortedY.end());
      subXLo = reInitializeSubsetYCor(&sortedX);
      subXHi = sortedX;

      treeNode* leftChild = kdtree_build_rec(subXLo, subYLo, treeDepth + 1);
      treeNode* rightChild = kdtree_build_rec(subXHi, subYHi, treeDepth + 1);

      treeNode *internalNode = initializeNode(internalPoint.x, internalPoint.y, depth, leftChild, rightChild);

      return internalNode;
      
    }

  }

}

/* treeHeight
 * 
 * Inputs:
 *     A treeNode pointer pointing to the root node in the tree.
 *
 * Return Values:
 *     This function returns the height of the tree being pointed to.
 *
 * Functional Description:
 *     This function recursively calls treeHeight on the two children of the root node, adding 1 to the child with the greatest height, ultimately
 *     propogating back up the tree, returning the max depth of the tree.
 */

int treeHeight(treeNode* node) {

  if (node == NULL) {
    return 0;
  }

  else {
    int lDepth = treeHeight(node->left);
    int rDepth = treeHeight(node->right);

    if (lDepth > rDepth) {
      return (lDepth + 1);
    }

    else {
      return (rDepth + 1);
    }

  }

}

/* kdtree_build
 * 
 * Inputs:
 *     A pointer to a vector of points.
 *
 * Return Values:
 *     This function returns whether or not the two points are the same.
 *
 * Functional Description:
 *     This function compares all the x and y cooridnates of the two points, and if all are the same, then the two segments are the same, 
 *     and the function returns true. This function returns false otherwise.
 */

kdtree* kdtree_build(vector<point2D> points) {

  vector<point2D> pointVector = points;

  if (pointVector.size() == 0) {
    return NULL;
  }

  kdtree* treePointer = new (nothrow) kdtree;
  assert(treePointer);
  vector<point2D> pointsByX, pointsByY;
  int treeDepth = 0;

  remove_coincident_points(&pointVector);

  pointsByX = pointVector;
  pointsByY = pointVector;

  sort(pointsByY.begin(), pointsByY.end(), pointSorterY);

  treePointer->root = kdtree_build_rec(pointsByX, pointsByY, treeDepth);
  treePointer->count = pointsByX.size();
  treePointer->height = treeHeight(treePointer->root);

  return treePointer;

}

