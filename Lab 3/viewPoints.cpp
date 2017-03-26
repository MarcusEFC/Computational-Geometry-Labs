#include "rtimer.h"
#include "kdtree.h"

#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <vector>
#include <new>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat cyan[3] = {0.0, 1.0, 1.0};

/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);
void reset(); 
void colorRegion(float x1, float y1, float x2, float y2);

/* global variables */
const int WINDOWSIZE = 500;
const int POINT_SIZE  = 6.0f;

//the array of n points
vector<point2D>  points;
int n;

// the kd-tree created with the points
kdtree* tree;

/* ****************************** */
/* initialize  the array of points stored in global variable points[] with random points */
void initialize_points_random() {
  
  int i;
  for (i=0; i<n; i++) {

    point2D newPoint;
    newPoint.x = (int)(.1*WINDOWSIZE)/2 + rand() % ((int)(.9*WINDOWSIZE));
    newPoint.y =  (int)(.1*WINDOWSIZE)/2 + rand() % ((int)(.9*WINDOWSIZE));

    points.push_back(newPoint);

  }
}


//This test case tests that out code can detect and remove coincident points.
void initialize_points_identical () {
  point2D newPoint;

  newPoint.x = 100;
  newPoint.y = 123;

  for (int i = 0; i < 4; i++) {
    points.push_back(newPoint);
  }

  newPoint.x = 50;
  newPoint.y = 75;

  for (int i = 0; i < 4; i++) {
    points.push_back(newPoint);
  }

  newPoint.x = 148;
  newPoint.y = 198;

  points.push_back(newPoint);

  newPoint.x = 100;
  newPoint.y = 129;

  points.push_back(newPoint);

  newPoint.x = 108;
  newPoint.y = 167;

  points.push_back(newPoint);

  newPoint.x = 234;
  newPoint.y = 198;

  points.push_back(newPoint);

  newPoint.x = 87;
  newPoint.y = 198;

  points.push_back(newPoint);

}

//This test case tests that out code works for the degenerative case that can result in infinite recursion.
void initialize_points_degenerative () {
  point2D newPoint;

  newPoint.x = 26;
  newPoint.y = 66;

  points.push_back(newPoint);

  newPoint.x = 36;
  newPoint.y = 66;

  points.push_back(newPoint);

  newPoint.x = 36;
  newPoint.y = 56;

  points.push_back(newPoint);

}

//This test case tests that out code works when all points have the same y coordinate.
void initialize_points_horizontal () {
  point2D newPoint;

  newPoint.x = WINDOWSIZE / 2;
  newPoint.y = WINDOWSIZE / 2;

  for (int i = 0; i < 4; i ++) {
    points.push_back(newPoint);
    newPoint.x += 50;
  }

}

//This test case tests that out code works when all points have the same x coordinate.
void initialize_points_vertical () {
  point2D newPoint;

  newPoint.x = WINDOWSIZE / 2;
  newPoint.y = WINDOWSIZE / 2;

  for (int i = 0; i < 4; i ++) {
    points.push_back(newPoint);
    newPoint.y += 50;
  }
  
}

//This test case was written by Bobby, Demi and Carolina.
void initialize_points_square(){
  point2D a,b,c,d;
  
  a.x=100;
  a.y=100;
  b.x=100;
  b.y=200;
  c.x=200;
  c.y=100;
  d.x=200;
  d.y=200;
  points.push_back(a);
  points.push_back(b);
  points.push_back(c);
  points.push_back(d);
}

//This test case was written by Drew.
void initialize_points_cross() {

  int i;

  point2D newPoint;

  for (i=0; i < n / 2; i++) {
    newPoint.y = WINDOWSIZE / 2;
    newPoint.x = rand() % (WINDOWSIZE / 2) + (WINDOWSIZE / 4);
    points.push_back(newPoint);
  }
  for (i = n / 2; i < n; i++) {
    newPoint.x = WINDOWSIZE / 2;
    newPoint.y = rand() % (WINDOWSIZE / 2) + (WINDOWSIZE / 4);
    points.push_back(newPoint);
  }
}

//This test case was written by Drew.
void initialize_points_right_angle() {
    
  point2D newPoint;
  
  newPoint.x = 2 * WINDOWSIZE / 3;
  newPoint.y = 2 * WINDOWSIZE / 3;

  points.push_back(newPoint);

  int i;

  for (i=1; i < n / 2; i++) {
    newPoint.y = 2 * WINDOWSIZE / 3;
    newPoint.x = rand() % (WINDOWSIZE / 3) + (WINDOWSIZE / 3);
    points.push_back(newPoint);
  }

  for (i = n / 2; i < n; i++) {
    newPoint.x = 2 * WINDOWSIZE / 3;
    newPoint.y = rand() % (WINDOWSIZE / 3) + (WINDOWSIZE / 3);
    points.push_back(newPoint);
  }

}

/* ****************************** */
/* print the array of points stored in global variable points[]*/
void print_points() {
    //assert(points);
    int i;
    printf("points: ");
    for (i=0; i<n; i++) {
        printf("[%d,%df] ", points[i].x, points[i].y);
    }
    printf("\n");
    fflush(stdout);  //flush stdout, weird sync happens when using gl thread
}

/* ****************************** */
int main(int argc, char** argv) {
    
    // read number of points from user
    if (argc!=2) {
        printf("usage: viewPoints <nbPoints>\n");
        exit(1);
    }
    
    n = atoi(argv[1]);
    printf("you entered n=%d\n", n);
    assert(n >0);

    initialize_points_random();
    
    srand(time(NULL));
    reset();
    
    /* initialize GLUT  */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
    glutInitWindowPosition(100,100);
    glutCreateWindow(argv[0]);
    
    /* register callback functions */
    glutDisplayFunc(display);
    glutKeyboardFunc(keypress);
    
    /* init GL */
    /* set background color black*/
    glClearColor(0, 0, 0, 0);
    
    /* circular points */
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glPointSize(POINT_SIZE);

    /* give control to event handler */
    glutMainLoop();
    return 0;
}




/* ****************************** */
/* draw a single point */
void draw_point(point2D point)
{
    glColor3fv(yellow);

    glBegin(GL_POINTS);
    glVertex2f(point.x, point.y);
    glEnd();
}


/* *****************ma************* */
/* draw a line between two points */
void draw_line(lineSegment2D line)
{
    glColor3fv(black);
    glLineWidth(1);
    
    glBegin(GL_LINES);
    glVertex2f(line.p1.x, line.p1.y);
    glVertex2f(line.p2.x, line.p2.y);
    glEnd();
}

/* ****************************** */
/* draw the array of points stored in global variable points[] 
   each point is drawn as a small square 
  
*/
void draw_points() {

  //const int R = 1;
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  //set color 
  glColor3fv(yellow);   
  
  //assert(points);
  int i;
  for (i=0; i<n; i++) {
    draw_point(points[i]); 
  }

}

/* ****************************** */
/* Colors in the region of the tree leaf.
 */

/* colorRegion
 * 
 * Inputs:
 *     One vertex of the rectangle and the opposite point of the rectangle to be colored.
 *
 * Return Values:
 *     N/A
 *
 * Functional Description:
 *     This function randomly picks a color, and colors the region specified by the two points passed as parameters.
 */

void colorRegion(float x1, float y1, float x2, float y2) {

  int colorNum = rand() % 4 + 1;

  switch (colorNum) {
    case 1: 
      glColor3fv(yellow);
      break;
    case 2: 
      glColor3fv(blue);
      break;
    case 3: 
      glColor3fv(red);
      break;
    case 4: 
      glColor3fv(white);
      break;
  }

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glRectf(x1, y1, x2, y2);
  glEnd();

}

/* ****************************** */
/* recursive draw function for drawing a tree rooted at the given node
 */
void draw_node(treeNode *node, float xBoundaryMin, float yBoundaryMin, float xBoundaryMax, float yBoundaryMax) {

  if (node == NULL) {
    return;
  }

  else if (node->type == 'l') { //Leaf
    colorRegion(xBoundaryMin, yBoundaryMin, xBoundaryMax, yBoundaryMax);
    return;

  }

  else if (node->type == 'v') { //Vertical Line Node

    float newXMinLeft = xBoundaryMin;
    float newXMaxLeft = node->p.x;

    float newXMinRight = node->p.x;
    float newXMaxRight = WINDOWSIZE;

    draw_node(node->left, newXMinLeft, yBoundaryMin, newXMaxLeft, yBoundaryMax);
    draw_node(node->right, newXMinRight, yBoundaryMin, newXMaxRight, yBoundaryMax);

  }

  else { //Horizontal Line Node

    float newYMinLeft = yBoundaryMin;
    float newYMaxLeft = node->p.y;

    float newYMinRight = node->p.y;
    float newYMaxRight = WINDOWSIZE;

    draw_node(node->left, xBoundaryMin, newYMinLeft, xBoundaryMax, newYMaxLeft);
    draw_node(node->right, xBoundaryMin, newYMinRight, xBoundaryMax, newYMaxRight);

  }

}

/* draw_lines
 * 
 * Inputs:
 *     - A pointer to a node in the kd tree.
 *     - The starting and ending point of the line to be drawn.
 *
 * Return Values:
 *     N/A
 *
 * Functional Description:
 *     This function draws all the lines in the kdtree recursively. This function first determines what kind of line the node is, draws the line, updates
 *     the new boundaries for the remaining nodes, and calls draw_lines on both the nodes left and right children.
 */

void draw_lines(treeNode *node, float xBoundaryMin, float yBoundaryMin, float xBoundaryMax, float yBoundaryMax) {

  if (node->type == 'v') { //Vertical Line Node

    lineSegment2D vLine;

    vLine.p1.x = node->p.x;
    vLine.p1.y = yBoundaryMin;
    vLine.p2.x = node->p.x;
    vLine.p2.y = yBoundaryMax;

    draw_line(vLine);

    float newXMinLeft = xBoundaryMin;
    float newXMaxLeft = node->p.x;

    float newXMinRight = node->p.x;
    float newXMaxRight = WINDOWSIZE;

    draw_lines(node->left, newXMinLeft, yBoundaryMin, newXMaxLeft, yBoundaryMax);
    draw_lines(node->right, newXMinRight, yBoundaryMin, newXMaxRight, yBoundaryMax);

  }

  else if (node->type == 'h') { //Horizontal Line Node

    lineSegment2D hLine;

    hLine.p1.x = xBoundaryMin;
    hLine.p1.y = node->p.y;
    hLine.p2.x = xBoundaryMax;
    hLine.p2.y = node->p.y;

    draw_line(hLine);

    float newYMinLeft = yBoundaryMin;
    float newYMaxLeft = node->p.y;

    float newYMinRight = node->p.y;
    float newYMaxRight = WINDOWSIZE;

    draw_lines(node->left, xBoundaryMin, newYMinLeft, xBoundaryMax, newYMaxLeft);
    draw_lines(node->right, xBoundaryMin, newYMinRight, xBoundaryMax, newYMaxRight);

  }

  else {
    return;
  }

}

/* ****************************** */
/* draw the kd-tree stored in the global variable kdTree
 */
void draw_kdtree()
{
    assert(tree);
    draw_node(tree->root, 0, 0, (float) WINDOWSIZE, (float) WINDOWSIZE);
}

/* ****************************** */
void display(void) {
    
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); //clear the matrix
    
    
    /* the default GL window is [-1,1]x[-1,1] with the origin in the
     center the points are in the range (0,0) to (WINSIZE,WINSIZE), so
     they need to be mapped to [-1,1]x [-1,1] */
    glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);
    glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0);

    //eventually we'll want to call the function that draws the kdtree
    //draw_kdtree();
   
    //for now we just draw the input points 
    draw_points();
    draw_node(tree->root, 0, 0, WINDOWSIZE, WINDOWSIZE);
    //draw_lines(tree->root, 0, 0, WINDOWSIZE, WINDOWSIZE);

 
    /* execute the drawing commands */
    glFlush();
}

void reset() {

  //re-initialize points
  points.clear(); 
  initialize_points_random();

  //free current tree
  if (tree != NULL) {
    kdtree_free(tree); 
  }
  
  Rtimer rt1;
  rt_start(rt1);

  tree = kdtree_build(points);
  rt_stop(rt1);
  char buf [1024];
  rt_sprint(buf,rt1);
  printf("time to generate kd-tree:  %s\n\n", buf);
  fflush(stdout);
  
  // print the tree
  kdtree_print(tree);

}

void testCaseReset (int caseNum) {

  points.clear();

  switch (caseNum) {
    case 1: 
      initialize_points_identical();
      break;
    case 2: 
      initialize_points_degenerative();
      break;
    case 3: 
      initialize_points_vertical();
      break;
    case 4: 
      initialize_points_horizontal();
      break;
    case 5:
      initialize_points_square();
      break;
    case 6:
      initialize_points_cross();
      break;
    case 7:
      initialize_points_right_angle();
      break;
  }

  if (tree != NULL) {
    kdtree_free(tree); 
  }

  Rtimer rt1;
  rt_start(rt1);
  tree = kdtree_build(points);
  rt_stop(rt1);
  char buf [1024];
  rt_sprint(buf,rt1);
  printf("time to generate kd-tree:  %s\n\n", buf);
  fflush(stdout);

  // print the tree
  kdtree_print(tree);

}

/* ****************************** */
void keypress(unsigned char key, int x, int y) {
    switch(key)
    {
        case ' ':
          reset();
          glutPostRedisplay();
          break;

        case '1':
          testCaseReset(1);
          glutPostRedisplay();
          break;

        case '2':
          testCaseReset(2);
          glutPostRedisplay();
          break;

        case '3':
          testCaseReset(3);
          glutPostRedisplay();
          break;

        case '4':
          testCaseReset(4);
          glutPostRedisplay();
          break;

        case '5':
          testCaseReset(5);
          glutPostRedisplay();
          break;

        case '6':
          testCaseReset(6);
          glutPostRedisplay();
          break;

        case '7':
          testCaseReset(7);
          glutPostRedisplay();
          break;

        case 'q':
          exit(0);
          break;
    }
}

/* Handler for window re-size event. Called back when the window first appears and
 whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    
    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);
    
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset
    gluOrtho2D(0.0, (GLdouble) width, 0.0, (GLdouble) height); 
}