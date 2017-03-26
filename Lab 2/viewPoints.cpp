/* view.c 

Marcus Christiansen

What it does:  

Draws a set of horizontal and vertical line segments in the default 2D
projection. Then computes their intersections using the line sweep
algorithm, and  simulates the algorithm as it runs.

*/

#include "geom.h"
#include "rtimer.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <vector> 
#include <map>

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

GLint fillmode = 0;

/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);
void timerfunc(); 

void initialize_segments_random(); 
void initialize_segments_horizontal();
void initialize_segments_random2();
void initialize_segments_overlap();
void initialize_segments_fourintersections();
void initialize_segments_overlap2(); 
void print_segments();
void print_points(); 

//renders the sweep line 
void draw_sweep_line();
//renders the active structure
void draw_active_structure();
//renders the intersection points 
void draw_intersection_points();

/* global variables */
const int WINDOWSIZE = 500; 

int init_case = 0; 
const int NB_TEST_CASES = 2; 

//NOTE: all the structures below need to be global so that they can be rendered

//current position of sweep line 
int sweep_line_x = 0; 

//number of segments requested by user  
int n; 

//the array of  segments
vector<segment2D>  segments;

//the intersections points of the segments 
vector<point2D> intpoints; 

//the events 
vector<sweepEvent> events; 

/* ************************************************** */
//fills global variable "segments" with n segments 
void initialize_segments() {

  switch (init_case)  {
      
    case 0: 
      initialize_segments_random(); 
      break;
      
    case 1: 
      initialize_segments_horizontal(); 
      break; 
      
    default: 
      initialize_segments_random(); 
    }

  init_case = (init_case+1) % NB_TEST_CASES;
  return; 
}

/* ************************************************** */
//fills global variable "segments" with n horizontal segments 
void initialize_segments_horizontal() {

  int i; 
  point2D a,b;
  segment2D s; 

  //clear the vector
  segments.clear(); 

  //a long horizontal segment 
  a.x = 1; 
  a.y = WINDOWSIZE/2; 
  b.x = WINDOWSIZE - 10; 
  b.y = a.y; 

  s.start = a; s.end = b; 
  segments.push_back(s);  

  //n-1 vertical segments 
  for (i=0; i<n-1; i++) {
    
    a.x = i*WINDOWSIZE/n; 
    a.y = WINDOWSIZE/2 - random() % ((int)(.4*WINDOWSIZE)); 
    b.x = a.x; 
    b.y = WINDOWSIZE/2 + random() % ((int)(.4*WINDOWSIZE)); 
    s.start = a; s.end = b; 
    segments.push_back(s); 
  }

}


/* ****************************** */
//fills global variable "segments" with n random segments 
void initialize_segments_random() {
  
  //clear the vector 
  segments.clear(); 

  int i; 
  point2D a, b; 
  segment2D s; 
  for (i=0; i<n; i++) {
    if (random()%2 == 0) {
      //horizontal segment
      a.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
      a.y =  (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
      b.y = a.y; 
      b.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 

      if (a.x < b.x) {
	s.start = a; s.end = b; 
      } else {
	s.start = b; s.end = a; 
      } 
 
    } else {
      //vertical segment 
      a.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
      b.x = a.x; 
      a.y = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
      b.y = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 

      if (a.y < b.y) {
	s.start = a; s.end = b; 
      } else {
	s.start = b; s.end = a; 
      }
    }

    //insert the segment in the array of segments 
    segments.push_back (s); 
  } //for i
}

//This test case tests that out code works for overlapping vertical segments.
void initialize_segments_overlap_vertical(){
  segments.clear();
  point2D a,b;
  segment2D s;

  for (int i = 1; i < 4; i++) {
    a.x = 0;
    b.x = WINDOWSIZE;
    a.y = 1 + i * WINDOWSIZE / 4;
    b.y = a.y;
    
    s.start = a;
    s.end = b;
    segments.push_back(s);
  }

  double k = 0.25;
    
  for (int i = 1; i < 4; i++) {
    for (int j = 1; j < 3; j++) {
      a.x = WINDOWSIZE*k;
      b.x = WINDOWSIZE*k;
      a.y = 0;
      b.y = WINDOWSIZE;
      
      s.start = a;
      s.end = b;
      segments.push_back(s);
    }

  k += 0.25;

  }
    
}

//This test case tests that out code works for overlapping horizontal segments.
void initialize_segments_overlap_horizontal(){
  segments.clear();
  point2D a,b;
  segment2D s;

  int k = 0;

  for (int i = 0; i < 10; i++) {
    a.x = 0;
    b.x = k + 50;
    a.y = WINDOWSIZE/2;
    b.y = a.y;
    
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 50;
  }
    
  for (int i = 1; i < 4; i++) {
    a.x = 1 + i * WINDOWSIZE / 4;
    b.x = a.x;
    a.y = 1;
    b.y = WINDOWSIZE - 1;
    
    s.start = a;
    s.end = b;
    segments.push_back(s);
  }
   
}

//This test case tests that out code can detect several intersections simultanesouly.
void initialize_segments_lattice(){
  segments.clear();
  point2D a,b;
  segment2D s;

  int k = 10;

  for (int i = 0; i < 49; i++) { //Horizontal
    a.x = 0;
    b.x = WINDOWSIZE;
    a.y = k;
    b.y = a.y;
        
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 10;

  }

  k = 10;

  for (int i = 0; i < 49; i++) { //Vertical
    a.x = k;
    b.x = k;
    a.y = 0;
    b.y = WINDOWSIZE;
        
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 10;

  }

    
}

//This test case tests that out code can remove the correct segment from the active structure when several segments have the same y-coordinate
void initialize_segments_rolling(){
  segments.clear();
  point2D a,b;
  segment2D s;

  int k = 0;

  for (int i = 0; i < 9; i++) { //Horizontal
    a.x = k;
    b.x = a.x + 100;
    a.y = WINDOWSIZE / 2;
    b.y = a.y;

    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 50;

  }

  k = 50;

  for (int i = 0; i < 9; i++) { //Vertical
    a.x = k;
    b.x = k;
    a.y = 0;
    b.y = WINDOWSIZE;
        
    s.start = a;
    s.end = b;
    segments.push_back(s);

    k += 50;

  }
 
}

//This case tests that our code can detect the correct intersections when several horizontal and vertical liens overlap at the same point.
//This test case was written by: James H. Lemkemeier.
void initialize_segments_fourintersections(){
    segments.clear();
    point2D a,b,c,d;
    a.x=200;
    a.y=150;
    b.x=250;
    b.y=150;
    segment2D temp;
    temp.start=a;
    temp.end=b;
    segments.push_back(temp);
    segments.push_back(temp);
    c.x=225;
    c.y=175;
    d.x=225;
    d.y=150;
    temp.start=d;
    temp.end=c;
    segments.push_back(temp);
    segments.push_back(temp);
}

//This case tests that our code can detect the correct amount of intersections with several overlaps.
//This test case was written by: Jeonguk Choi.
void initialize_segments_overlap() {
    segments.clear();
    int i;
    point2D a,b;
    segment2D s;
    
    for (i = 0; i < 10; i++) {
        a.x = 1 + i * (WINDOWSIZE/20);
        b.x = WINDOWSIZE - 1 - i*(WINDOWSIZE/20);
        a.y = WINDOWSIZE/2;
        b.y = a.y;
        
        s.start = a;
        s.end = b;
        segments.push_back(s);
    }
    
    for (i = 1; i < 4; i++) {
        a.x = 1 + i * WINDOWSIZE / 4;
        b.x = a.x;
        a.y = 1;
        b.y = WINDOWSIZE - 1;
        
        s.start = a;
        s.end = b;
        segments.push_back(s);
    }
}

/* ************************************************** */
void print_segments() {

  for (int i=0; i<segments.size(); i++) {
    printf("segment %d: [(%d,%d), (%d,%d)]\n",
	   i, segments[i].start.x, segments[i].start.y, segments[i].end.x, segments[i].end.y);

  }

}

void print_points() {

  for (int i = 0; i < intpoints.size(); i++) {
    printf("Intersection %d: [X: %d, Y: %d]\n",
     i, intpoints[i].x, intpoints[i].y);

  }

}

/* ****************************** */
int main(int argc, char** argv) {

  //read number of points from user
  if (argc!=2) {
    printf("usage: viewPoints <nbPoints>\n");
    exit(1); 
  }
  n = atoi(argv[1]); 
  printf("you entered n=%d\n", n);
  assert(n >0); 

  initialize_segments_random();
  print_segments();

  Rtimer rt1; 
  rt_start(rt1); 

  intpoints = line_segment_scan(segments, events);

  rt_stop(rt1); 

  print_points();

  //print the timing 
  char buf [1024]; 
  rt_sprint(buf,rt1);
  printf("run time:  %s\n\n", buf);
  fflush(stdout); 

  /* initialize GLUT  */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display); 
  glutKeyboardFunc(keypress);
  glutIdleFunc(timerfunc); 

  /* init GL */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);   
  /* here we can enable depth testing and double buffering and so
     on */

  /* give control to event handler */
  glutMainLoop();
  return 0;
}


/* ****************************** */
/* draw the segments stored in global variable segments */
void draw_segments(){

  //set color 
  glColor3fv(yellow);   
  
  int i;
  for (i=0; i<segments.size(); i++) {
    glBegin(GL_LINES);
    glVertex2f(segments[i].start.x, segments[i].start.y); 
    glVertex2f(segments[i].end.x, segments[i].end.y);
    glEnd();
  }

}

//draw the sweep line 
void draw_sweep_line() {

  //sweep line color 
  glColor3fv(green); 

  //the current position of sweep line is sweep_line_x; assume it's a
  //segment from y=0 to y=windowsize;
  glBegin(GL_LINES); 
  glVertex2f(sweep_line_x, 0); 
  glVertex2f(sweep_line_x, WINDOWSIZE); 
  glEnd();
}

//draw all the elements in the active structure 
void draw_active_structure() {

  for (int i = 0; i < segments.size(); i++) {
    if (segments[i].end.x <= sweep_line_x) {

      glColor3fv(yellow);
      glBegin(GL_LINES);
      glVertex2f(segments[i].start.x, segments[i].start.y); 
      glVertex2f(segments[i].end.x, segments[i].end.y);
      glEnd();

    }

    else if (segments[i].start.x <= sweep_line_x) {

      glColor3fv(red);
      glBegin(GL_LINES);
      glVertex2f(segments[i].start.x, segments[i].start.y); 
      glVertex2f(segments[i].end.x, segments[i].end.y);
      glEnd();

    }

  }
  
}

//draw all the elements in intpoints 
void draw_intersection_points() {

  for (int i = 0; i < intpoints.size(); i++) {
    if (intpoints[i].x <= sweep_line_x) {

      glEnable( GL_POINT_SMOOTH );
      glColor3fv(blue);  
      glPointSize(3.5);
      glBegin(GL_POINTS);
      glVertex2f(intpoints[i].x, intpoints[i].y);
      glEnd();

    }

  }

}

/* ****************************** */
void display(void) {

  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); //clear the matrix


  /* The default GL window is [-1,1]x[-1,1] with the origin in the
     center. 
     
     The points are in the range (0,0) to (WINSIZE,WINSIZE), so they
     need to be mapped to [-1,1]x [-1,1]x */
  glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);  
  glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0); 

  draw_segments();
  draw_active_structure(); 
  draw_intersection_points(); 
  draw_sweep_line(); 

  /* execute the drawing commands */
  glFlush();
}

/* ****************************** */
void keypress(unsigned char key, int x, int y) {

  Rtimer rt1;
  char buf [1024];

  switch(key) {
  case 'q':
    exit(0);
    break;

  case 'i': 
    initialize_segments(); 
    glutPostRedisplay();
    sweep_line_x = 0;

    rt_start(rt1); 
    intpoints = line_segment_scan(segments, events);
    rt_stop(rt1); 

    print_segments();
    print_points();

    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout); 
    break;

  case '1': 
    initialize_segments_overlap_vertical(); 
    glutPostRedisplay();
    sweep_line_x = 0;

    rt_start(rt1); 
    intpoints = line_segment_scan(segments, events);
    rt_stop(rt1); 

    print_segments();
    print_points();

    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout);  
    break;

  case '2': 
    initialize_segments_overlap_horizontal(); 
    glutPostRedisplay();
    sweep_line_x = 0;

    rt_start(rt1); 
    intpoints = line_segment_scan(segments, events);
    rt_stop(rt1);

    print_segments();
    print_points(); 

    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout);  
    break;

  case '3': 
    initialize_segments_lattice(); 
    glutPostRedisplay();
    sweep_line_x = 0;

    rt_start(rt1); 
    intpoints = line_segment_scan(segments, events);
    rt_stop(rt1);

    print_segments();
    print_points(); 

    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout);  
    break;

  case '4': 
    initialize_segments_rolling(); 
    glutPostRedisplay();
    sweep_line_x = 0;

    rt_start(rt1); 
    intpoints = line_segment_scan(segments, events);
    rt_stop(rt1);

    print_segments();
    print_points(); 

    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout);  
    break;

  case '5': 
    initialize_segments_fourintersections(); 
    glutPostRedisplay();
    sweep_line_x = 0;

    rt_start(rt1); 
    intpoints = line_segment_scan(segments, events);
    rt_stop(rt1); 

    print_segments();
    print_points();

    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout);   
    break;

  case '6': 
    initialize_segments_overlap(); 
    glutPostRedisplay();
    sweep_line_x = 0;

    rt_start(rt1); 
    intpoints = line_segment_scan(segments, events);
    rt_stop(rt1); 

    print_segments();
    print_points();

    rt_sprint(buf,rt1);
    printf("run time:  %s\n\n", buf);
    fflush(stdout);  
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

void timerfunc() {
  
  /* LT: I used this to slow things down in a controlled way; probably not
     necessary for this assignment..*/
  static int lastFrameTime=0;  
  //note: a static variable, remembered from one call to the next
  int now, elapsed_ms; 
  
  now = glutGet (GLUT_ELAPSED_TIME); 
  elapsed_ms = now - lastFrameTime; 
  lastFrameTime=now; 
  
  sweep_line_x++; 

  
  glutPostRedisplay(); 

}
