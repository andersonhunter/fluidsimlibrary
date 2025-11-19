// openGL and Glut adapted from CS450 starter code from Professor Mike Bailey, Oregon State University

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <chrono>

#ifndef F_PI
#define F_PI		((float)(M_PI))
#define F_2_PI		((float)(2.f*F_PI))
#define F_PI_2		((float)(F_PI/2.f))
#endif


#ifdef WIN32
#include <windows.h>
#pragma warning(disable:4996)
#endif


// Decide if we're using the minmod limiter
#ifndef MINMOD
#define MINMOD false
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include "glew.h"
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#define MSEC	10000

#ifndef SIZE
#define SIZE	16
#endif

#ifndef CELLSIZE
#define CELLSIZE (1.0 / (float)SIZE)
#endif

// Number of Jacobi iterations to perform
// More iterations = more resolution (but more computation time ):
#ifndef JACOBIS
#define JACOBIS 1e-6
#endif

#include "glut.h"

#ifndef DEBUG
#define DEBUG true
#endif

struct Point {
	// A point within the grid which holds different quantities
	float temperature;  // Temperature
	float vx;           // X-velocity
	float vy;           // Y-Velocity
	float oldPressure;  // Previous pressure for Poisson
	float pressure;     // Temporary pressure, should be ~0 after projection
	float density;      // Density
	float sourceTerm;   // Divergence at the given point
};

uint64_t currentTime;

//	This is a sample OpenGL / GLUT program
//
//	The objective is to draw a 3d object and change the color of the axes
//		with a glut menu
//
//	The left mouse button does rotation
//	The middle mouse button does scaling
//	The user interface allows:
//		1. The axes to be turned on and off
//		2. The color of the axes to be changed
//		3. Debugging to be turned on and off
//		4. Depth cueing to be turned on and off
//		5. The projection to be changed
//		6. The transformations to be reset
//		7. The program to quit
//
//	Author:			Joe Graphics

// title of these windows:

const char *WINDOWTITLE = "Fluid Simulation -- Hunter Anderson";
const char *GLUITITLE   = "User Interface Window";

// what the glui package defines as true and false:

const int GLUITRUE  = true;
const int GLUIFALSE = false;

// the escape key:

const int ESCAPE = 0x1b;

// initial window size:

const int INIT_WINDOW_SIZE = 1000;

// size of the 3d box to be drawn:

const float BOXSIZE = 2.f;

// Size of the grid
const float GRIDCELLSIZE = 2.;

// Proportion of the arrow size that the lil points comprise
const float ARROWFRAC = 0.6;
const float ARROWOFF = 0.2; // Offset of arrow within cell

// multiplication factors for input interaction:
//  (these are known from previous experience)

const float ANGFACT = 1.f;
const float SCLFACT = 0.005f;

// minimum allowable scale factor:

const float MINSCALE = 0.05f;

// scroll wheel button values:

const int SCROLL_WHEEL_UP   = 3;
const int SCROLL_WHEEL_DOWN = 4;

// equivalent mouse movement when we click the scroll wheel:

const float SCROLL_WHEEL_CLICK_FACTOR = 5.f;

// active mouse buttons (or them together):

const int LEFT   = 4;
const int MIDDLE = 2;
const int RIGHT  = 1;

// which projection:

enum Projections
{
	ORTHO,
	PERSP
};

// which button:

enum ButtonVals
{
	RESET,
	QUIT
};

// window background color (rgba):

const GLfloat BACKCOLOR[ ] = { 0., 0., 0., 1. };

// line width for the axes:

const GLfloat AXES_WIDTH   = 3.;

// the color numbers:
// this order must match the radio button order, which must match the order of the color names,
// 	which must match the order of the color RGB values

enum Colors
{
	RED,
	YELLOW,
	GREEN,
	CYAN,
	BLUE,
	MAGENTA
};

char * ColorNames[ ] =
{
	(char *)"Red",
	(char*)"Yellow",
	(char*)"Green",
	(char*)"Cyan",
	(char*)"Blue",
	(char*)"Magenta"
};

// the color definitions:
// this order must match the menu order

const GLfloat Colors[ ][3] = 
{
	{ 1., 0., 0. },		// red
	{ 1., 1., 0. },		// yellow
	{ 0., 1., 0. },		// green
	{ 0., 1., 1. },		// cyan
	{ 0., 0., 1. },		// blue
	{ 1., 0., 1. },		// magenta
};

// fog parameters:

const GLfloat FOGCOLOR[4] = { .0f, .0f, .0f, 1.f };
const GLenum  FOGMODE     = GL_LINEAR;
const GLfloat FOGDENSITY  = 0.30f;
const GLfloat FOGSTART    = 1.5f;
const GLfloat FOGEND      = 4.f;

// for lighting:

const float	WHITE[ ] = { 1.,1.,1.,1. };

// for animation:

const int MS_PER_CYCLE = 10000;		// 10000 milliseconds = 10 seconds

// non-constant global variables:

int		ActiveButton;			// current button that is down
GLuint	AxesList;				// list to hold the axes
int		AxesOn;					// != 0 means to draw the axes
GLuint	GridList;				// object display list
GLuint	ArrowList;
GLuint	BackgroundList;
int		DebugOn;				// != 0 means to print debugging info
int		DepthCueOn;				// != 0 means to use intensity depth cueing
int		DepthBufferOn;			// != 0 means to use the z-buffer
int		DepthFightingOn;		// != 0 means to force the creation of z-fighting
int		MainWindow;				// window id for main graphics window
int		NowColor;				// index into Colors[ ]
int		NowProjection;		// ORTHO or PERSP
float	Scale;					// scaling factor
int		ShadowsOn;				// != 0 means to turn shadows on
float	Time;					// used for animation, this has a value between 0. and 1.
int		Xmouse, Ymouse;			// mouse values
float	Xrot, Yrot;				// rotation angles in degrees
struct Point* grid;
float	maxMagnitude = 0.000001;

// function prototypes:

void	Animate( );
void	Display( );
void	DoAxesMenu( int );
void	DoColorMenu( int );
void	DoDepthBufferMenu( int );
void	DoDepthFightingMenu( int );
void	DoDepthMenu( int );
void	DoDebugMenu( int );
void	DoMainMenu( int );
void	DoProjectMenu( int );
void	DoRasterString( float, float, float, char * );
void	DoStrokeString( float, float, float, float, char * );
float	ElapsedSeconds( );
void	InitGraphics( );
void	InitLists( );
void	InitMenus( );
void	Keyboard( unsigned char, int, int );
void	MouseButton( int, int, int, int );
void	MouseMotion( int, int );
void	Reset( );
void	Resize( int, int );
void	Visibility( int );

uint64_t timeSinceEpochMillisec();

void			Axes( float );
void			HsvRgb( float[3], float [3] );
void			Cross(float[3], float[3], float[3]);
float			Dot(float [3], float [3]);
float			Unit(float [3], float [3]);
float			Unit(float [3]);

// Function prototypes
void calculateAdvection(float timestep);
struct Point getAtIndex(int x, int y);
void setVxAtIndex(int x, int y, float vx);
void setVyAtIndex(int x, int y, float vy);
float cubicInterpolate(float p0, float p1, float p2, float p3, float dt);
void InitGrid();
float solveDivergence(float p0, float p1, float p2, float p3, float sourceTerm);
float getSafePressure(float* pressureGrid, int x, int y);
float minmod(float a, float b);
void InitBoundaries();

// utility to create an array from 3 separate values:

float *
Array3( float a, float b, float c )
{
	static float array[4];

	array[0] = a;
	array[1] = b;
	array[2] = c;
	array[3] = 1.;
	return array;
}

// utility to create an array from a multiplier and an array:

float *
MulArray3( float factor, float array0[ ] )
{
	static float array[4];

	array[0] = factor * array0[0];
	array[1] = factor * array0[1];
	array[2] = factor * array0[2];
	array[3] = 1.;
	return array;
}


float *
MulArray3(float factor, float a, float b, float c )
{
	static float array[4];

	float* abc = Array3(a, b, c);
	array[0] = factor * abc[0];
	array[1] = factor * abc[1];
	array[2] = factor * abc[2];
	array[3] = 1.;
	return array;
}


float
Ranf( float low, float high )
{
        float r = (float) rand();               // 0 - RAND_MAX
        float t = r  /  (float) RAND_MAX;       // 0. - 1.

        return   low  +  t * ( high - low );
}

// call this if you want to force your program to use
// a different random number sequence every time you run it:
void
TimeOfDaySeed( )
{
	struct tm y2k;
	y2k.tm_hour = 0;    y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 2000; y2k.tm_mon = 0; y2k.tm_mday = 1;

	time_t  now;
	time( &now );
	double seconds = difftime( now, mktime(&y2k) );
	unsigned int seed = (unsigned int)( 1000.*seconds );    // milliseconds
	srand( seed );
}

// these are here for when you need them -- just uncomment the ones you need:

#include "setmaterial.cpp"
#include "setlight.cpp"
//#include "osusphere.cpp"
//#include "osucube.cpp"
//#include "osucylindercone.cpp"
//#include "osutorus.cpp"
//#include "bmptotexture.cpp"
//#include "loadobjmtlfiles.cpp"
#include "keytime.cpp"
//#include "glslprogram.cpp"
//#include "vertexbufferobject.cpp"
Keytimes Xpos;

// main program:

int
main( int argc, char *argv[ ] )
{
	// turn on the glut package:
	// (do this before checking argc and argv since glutInit might
	// pull some command line arguments out)

	glutInit( &argc, argv );

	// setup all the graphics stuff:

	InitGraphics( );

	// create the display lists that **will not change**:

	InitLists( );

	// init all the global variables used by Display( ):
	// this will also post a redisplay

	Reset( );

	// setup all the user interface stuff:

	InitMenus( );

	// Set up grid
	InitGrid();
	// Initialize boundaries
	InitBoundaries();

	currentTime = timeSinceEpochMillisec();

	// draw the scene once and wait for some interaction:
	// (this will never return)

	glutSetWindow( MainWindow );
	glutMainLoop( );

	// glutMainLoop( ) never actually returns
	// the following line is here to make the compiler happy:

	return 0;
}


// this is where one would put code that is to be called
// everytime the glut main loop has nothing to do
//
// this is typically where animation parameters are set
//
// do not call Display( ) from here -- let glutPostRedisplay( ) do it

void
Animate( )
{
	// put animation stuff in here -- change some global variables for Display( ) to find:

	int ms = glutGet(GLUT_ELAPSED_TIME);
	ms %= MS_PER_CYCLE;							// makes the value of ms between 0 and MS_PER_CYCLE-1
	Time = (float)ms / (float)MS_PER_CYCLE;		// makes the value of Time between 0. and slightly less than 1.

	// for example, if you wanted to spin an object in Display( ), you might call: glRotatef( 360.f*Time,   0., 1., 0. );

	// force a call to Display( ) next time it is convenient:
	// Compute max velocity
	float maxVelocity = 0.;
	for (int i = 0; i < SIZE * SIZE; i++) {
		float vx = grid[i].vx;
		float vy = grid[i].vy;
		float magnitude = sqrtf(vx * vx + vy * vy);
		if (magnitude > maxVelocity) {
			maxVelocity = magnitude;
		}
	}

	// CFL compliance
	float dt;
	const float Cmax = 0.5;
	const float kinematicViscosity = 0.001;

	if (maxVelocity > 0.f) {
		// CFL limit for advection
		float dtCFL = Cmax * ((float)CELLSIZE / maxVelocity);

		// Set explicit diffusion stability limit
		float dtVisc = 0.25 * ((float)CELLSIZE * (float)CELLSIZE) / kinematicViscosity;
		
		// Choose the smaller timestep
		dt = fmin(dtCFL, dtVisc);
	}
	else {
		// Choose reasonable timestep if grid is still
		dt = 0.01;
	}

	maxMagnitude = (maxVelocity >= 0.f) ? maxVelocity : 0.0001f;

	fprintf(stderr, "dt = %.6f, maxVelocity = %.5f\n", dt, maxVelocity);

	calculateAdvection(dt);
	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// draw the complete scene:

void
Display( )
{
	if (DebugOn != 0)
		fprintf(stderr, "Starting Display.\n");

	// set which window we want to do the graphics into:
	glutSetWindow( MainWindow );

	// erase the background:
	glDrawBuffer( GL_BACK );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable( GL_DEPTH_TEST );

	// specify shading to be flat:

	glShadeModel( GL_SMOOTH );

	// set the viewport to be a square centered in the window:

	GLsizei vx = glutGet( GLUT_WINDOW_WIDTH );
	GLsizei vy = glutGet( GLUT_WINDOW_HEIGHT );
	GLsizei v = vx < vy ? vx : vy;			// minimum dimension
	GLint xl = ( vx - v ) / 2;
	GLint yb = ( vy - v ) / 2;
	glViewport( xl, yb,  v, v );


	// set the viewing volume:
	// remember that the Z clipping  values are given as DISTANCES IN FRONT OF THE EYE
	// USE gluOrtho2D( ) IF YOU ARE DOING 2D !

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity( );
	if( NowProjection == ORTHO )
		glOrtho(-25.f, 25.f, -25.f, 25.f, 0.1f, 1000.f);
	else
		gluPerspective( 70.f, 1.f,	0.1f, 1000.f );

	// place the objects into the scene:

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );

	// set the eye position, look-at position, and up-vector:

	gluLookAt( 0.f, 30.f, 3.f,     0.f, 0.f, 0.f,     0.f, 1.f, 0.f );

	// rotate the scene:

	glRotatef( (GLfloat)Yrot, 0.f, 1.f, 0.f );
	glRotatef( (GLfloat)Xrot, 1.f, 0.f, 0.f );

	// uniformly scale the scene:

	if( Scale < MINSCALE )
		Scale = MINSCALE;
	glScalef( (GLfloat)Scale, (GLfloat)Scale, (GLfloat)Scale );

	// set the fog parameters:

	if( DepthCueOn != 0 )
	{
		glFogi( GL_FOG_MODE, FOGMODE );
		glFogfv( GL_FOG_COLOR, FOGCOLOR );
		glFogf( GL_FOG_DENSITY, FOGDENSITY );
		glFogf( GL_FOG_START, FOGSTART );
		glFogf( GL_FOG_END, FOGEND );
		glEnable( GL_FOG );
	}
	else
	{
		glDisable( GL_FOG );
	}

	// possibly draw the axes:

	if( AxesOn != 0 )
	{
		glColor3fv( &Colors[NowColor][0] );
		glCallList( AxesList );
	}

	// since we are using glScalef( ), be sure the normals get unitized:

	glEnable( GL_NORMALIZE );


	// draw the box object by calling up its display list:
	SetPointLight(GL_LIGHT0, 0., 1., 0., 1., 1., 1.);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);

	int msec = glutGet(GLUT_ELAPSED_TIME) % MSEC;

	float nowTime = (float)msec / 1000.f;

	glPushMatrix();
	glCallList(GridList);
	glPopMatrix();

	// Rotate, then transform
	for (int row = 0; row < SIZE; row++) {
		for (int col = 0; col < SIZE; col++) {
			glPushMatrix();
			glTranslatef(GRIDCELLSIZE * (0.5 - 0.5 * (float)SIZE + col), 0., GRIDCELLSIZE * (0.5 - 0.5 * (float)SIZE + row));
			struct Point nowPoint = getAtIndex(row, col);
			float vx = nowPoint.vx;
			float vy = nowPoint.vy;
			float theta = atan2(vx, vy);
			glRotatef(theta * 180. / M_PI, 0., 1., 0.);
			// Scale the arrow to [0, 1] based on magnitude
			float mag = sqrtf((vx * vx) + (vy * vy)) / maxMagnitude;
			glScalef(mag, mag, mag);
			glCallList(ArrowList);
			glPopMatrix();
		}
	}

	glPushMatrix();
	glCallList(BackgroundList);
	glPopMatrix();

	// swap the double-buffered framebuffers:

	glutSwapBuffers( );

	// be sure the graphics buffer has been sent:
	// note: be sure to use glFlush( ) here, not glFinish( ) !

	glFlush( );
}


void
DoAxesMenu( int id )
{
	AxesOn = id;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


void
DoColorMenu( int id )
{
	NowColor = id - RED;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


void
DoDebugMenu( int id )
{
	DebugOn = id;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


void
DoDepthBufferMenu( int id )
{
	DepthBufferOn = id;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


void
DoDepthFightingMenu( int id )
{
	DepthFightingOn = id;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


void
DoDepthMenu( int id )
{
	DepthCueOn = id;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// main menu callback:

void
DoMainMenu( int id )
{
	switch( id )
	{
		case RESET:
			Reset( );
			break;

		case QUIT:
			// gracefully close out the graphics:
			// gracefully close the graphics window:
			// gracefully exit the program:
			glutSetWindow( MainWindow );
			glFinish( );
			glutDestroyWindow( MainWindow );
			exit( 0 );
			break;

		default:
			fprintf( stderr, "Don't know what to do with Main Menu ID %d\n", id );
	}

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


void
DoProjectMenu( int id )
{
	NowProjection = id;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// use glut to display a string of characters using a raster font:

void
DoRasterString( float x, float y, float z, char *s )
{
	glRasterPos3f( (GLfloat)x, (GLfloat)y, (GLfloat)z );

	char c;			// one character to print
	for( ; ( c = *s ) != '\0'; s++ )
	{
		glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_24, c );
	}
}


// use glut to display a string of characters using a stroke font:

void
DoStrokeString( float x, float y, float z, float ht, char *s )
{
	glPushMatrix( );
		glTranslatef( (GLfloat)x, (GLfloat)y, (GLfloat)z );
		float sf = ht / ( 119.05f + 33.33f );
		glScalef( (GLfloat)sf, (GLfloat)sf, (GLfloat)sf );
		char c;			// one character to print
		for( ; ( c = *s ) != '\0'; s++ )
		{
			glutStrokeCharacter( GLUT_STROKE_ROMAN, c );
		}
	glPopMatrix( );
}


// return the number of seconds since the start of the program:

float
ElapsedSeconds( )
{
	// get # of milliseconds since the start of the program:

	int ms = glutGet( GLUT_ELAPSED_TIME );

	// convert it to seconds:

	return (float)ms / 1000.f;
}



// initialize the glut and OpenGL libraries:
//	also setup callback functions

void
InitGraphics( )
{
	if (DebugOn != 0)
		fprintf(stderr, "Starting InitGraphics.\n");

	Xpos.Init();
	Xpos.AddTimeValue(0.0, 0.000);
	Xpos.AddTimeValue(0.5, 2.718);
	Xpos.AddTimeValue(2.0, 0.333);
	Xpos.AddTimeValue(5.0, 3.142);
	Xpos.AddTimeValue(8.0, 2.718);
	Xpos.AddTimeValue(10.0, 0.000);

	// request the display modes:
	// ask for red-green-blue-alpha color, double-buffering, and z-buffering:

	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );

	// set the initial window configuration:

	glutInitWindowPosition( 0, 0 );
	glutInitWindowSize( INIT_WINDOW_SIZE, INIT_WINDOW_SIZE );

	// open the window and set its title:

	MainWindow = glutCreateWindow( WINDOWTITLE );
	glutSetWindowTitle( WINDOWTITLE );

	// set the framebuffer clear values:

	glClearColor( BACKCOLOR[0], BACKCOLOR[1], BACKCOLOR[2], BACKCOLOR[3] );

	// setup the callback functions:
	// DisplayFunc -- redraw the window
	// ReshapeFunc -- handle the user resizing the window
	// KeyboardFunc -- handle a keyboard input
	// MouseFunc -- handle the mouse button going down or up
	// MotionFunc -- handle the mouse moving with a button down
	// PassiveMotionFunc -- handle the mouse moving with a button up
	// VisibilityFunc -- handle a change in window visibility
	// EntryFunc	-- handle the cursor entering or leaving the window
	// SpecialFunc -- handle special keys on the keyboard
	// SpaceballMotionFunc -- handle spaceball translation
	// SpaceballRotateFunc -- handle spaceball rotation
	// SpaceballButtonFunc -- handle spaceball button hits
	// ButtonBoxFunc -- handle button box hits
	// DialsFunc -- handle dial rotations
	// TabletMotionFunc -- handle digitizing tablet motion
	// TabletButtonFunc -- handle digitizing tablet button hits
	// MenuStateFunc -- declare when a pop-up menu is in use
	// TimerFunc -- trigger something to happen a certain time from now
	// IdleFunc -- what to do when nothing else is going on

	glutSetWindow( MainWindow );
	glutDisplayFunc( Display );
	glutReshapeFunc( Resize );
	glutKeyboardFunc( Keyboard );
	glutMouseFunc( MouseButton );
	glutMotionFunc( MouseMotion );
	glutPassiveMotionFunc(MouseMotion);
	//glutPassiveMotionFunc( NULL );
	glutVisibilityFunc( Visibility );
	glutEntryFunc( NULL );
	glutSpecialFunc( NULL );
	glutSpaceballMotionFunc( NULL );
	glutSpaceballRotateFunc( NULL );
	glutSpaceballButtonFunc( NULL );
	glutButtonBoxFunc( NULL );
	glutDialsFunc( NULL );
	glutTabletMotionFunc( NULL );
	glutTabletButtonFunc( NULL );
	glutMenuStateFunc( NULL );
	glutTimerFunc( -1, NULL, 0 );

	// setup glut to call Animate( ) every time it has
	// 	nothing it needs to respond to (which is most of the time)
	// we don't need to do this for this program, and really should set the argument to NULL
	// but, this sets us up nicely for doing animation

	glutIdleFunc( Animate );

	// init the glew package (a window must be open to do this):

#ifdef WIN32
	GLenum err = glewInit( );
	if( err != GLEW_OK )
	{
		fprintf( stderr, "glewInit Error\n" );
	}
	else
		fprintf( stderr, "GLEW initialized OK\n" );
	fprintf( stderr, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
#endif

	// all other setups go here, such as GLSLProgram and KeyTime setups:

}


// initialize the display lists that will not change:
// (a display list is a way to store opengl commands in
//  memory so that they can be played back efficiently at a later time
//  with a call to glCallList( )

void
InitLists( )
{
	if (DebugOn != 0)
		fprintf(stderr, "Starting InitLists.\n");

	// ListName = glGenLists(1);
	// glNewList(ListName, GL_COMPILE);
	// glBegin(SHAPE);
	// glVertex3f(x, y, z);
	// glEnd();
	// glEndList();

	float xPoint = GRIDCELLSIZE * ((float)SIZE / 2);

	// Create the grid
	GridList = glGenLists(1);
	glNewList(GridList, GL_COMPILE);
		glBegin(GL_LINES);
		// Draw horizontal gridlines
		for (int i = 0; i <= SIZE; i++) {
			glColor3f(1., 0., 0.);
			// Draw the horizontal line
			glVertex3f(-xPoint, 0., xPoint - (GRIDCELLSIZE * (float)i));
			glVertex3f(xPoint, 0., xPoint - (GRIDCELLSIZE * (float)i));
			// Draw the vertical line
			glVertex3f(xPoint - (GRIDCELLSIZE * (float)i), 0., -xPoint);
			glVertex3f(xPoint - (GRIDCELLSIZE * (float)i), 0., xPoint);
		}
		glEnd();
	glEndList();

	// Draw the arrows
	float cf2 = GRIDCELLSIZE / 2.;
	ArrowList = glGenLists(1);
	glNewList(ArrowList, GL_COMPILE);
	glBegin(GL_LINE_STRIP);
	glColor3f(0., 0., 1.);
	// Left side of arrow
	glVertex3f(-cf2 + ARROWOFF, 0., 0.);
	// R side of arrow
	glVertex3f(cf2 - ARROWOFF, 0., 0.);
	// Bottom point of arrow
	glVertex3f(cf2 - (ARROWOFF + ARROWFRAC), 0., ARROWFRAC);
	// Back to the right side
	glVertex3f(cf2 - ARROWOFF, 0., 0.);
	// Aaaaand up to the top
	glVertex3f(cf2 - (ARROWFRAC + ARROWOFF), 0., -ARROWFRAC);
	glEnd();
	glEndList();

	// Draw a backboard for the grid
	BackgroundList = glGenLists(1);
	glNewList(BackgroundList, GL_COMPILE);
		glBegin(GL_QUADS);
		glColor3f(0.8, 0.8, 0.8);
		glVertex3f(-xPoint - GRIDCELLSIZE, -0.5, -xPoint - (2. * GRIDCELLSIZE));
		glVertex3f(xPoint + GRIDCELLSIZE, -0.5, -xPoint - (2. * GRIDCELLSIZE));
		glVertex3f(xPoint + GRIDCELLSIZE, -0.5, xPoint + 2. * GRIDCELLSIZE);
		glVertex3f(-xPoint - GRIDCELLSIZE, -0.5, xPoint + 2. * GRIDCELLSIZE);
		glEnd();
	glEndList();

	AxesList = glGenLists( 1 );
	glNewList( AxesList, GL_COMPILE );
		glLineWidth( AXES_WIDTH );
			Axes( 1.5 );
		glLineWidth( 1. );
	glEndList( );
}


// initialize the glui window:

void
InitMenus( )
{
	if (DebugOn != 0)
		fprintf(stderr, "Starting InitMenus.\n");

	glutSetWindow( MainWindow );

	int numColors = sizeof( Colors ) / ( 3*sizeof(float) );
	int colormenu = glutCreateMenu( DoColorMenu );
	for( int i = 0; i < numColors; i++ )
	{
		glutAddMenuEntry( ColorNames[i], i );
	}

	int axesmenu = glutCreateMenu( DoAxesMenu );
	glutAddMenuEntry( "Off",  0 );
	glutAddMenuEntry( "On",   1 );

	int depthcuemenu = glutCreateMenu( DoDepthMenu );
	glutAddMenuEntry( "Off",  0 );
	glutAddMenuEntry( "On",   1 );

	int depthbuffermenu = glutCreateMenu( DoDepthBufferMenu );
	glutAddMenuEntry( "Off",  0 );
	glutAddMenuEntry( "On",   1 );

	int depthfightingmenu = glutCreateMenu( DoDepthFightingMenu );
	glutAddMenuEntry( "Off",  0 );
	glutAddMenuEntry( "On",   1 );

	int debugmenu = glutCreateMenu( DoDebugMenu );
	glutAddMenuEntry( "Off",  0 );
	glutAddMenuEntry( "On",   1 );

	int projmenu = glutCreateMenu( DoProjectMenu );
	glutAddMenuEntry( "Orthographic",  ORTHO );
	glutAddMenuEntry( "Perspective",   PERSP );

	int mainmenu = glutCreateMenu( DoMainMenu );
	glutAddSubMenu(   "Axes",          axesmenu);
	glutAddSubMenu(   "Axis Colors",   colormenu);

#ifdef DEMO_DEPTH_BUFFER
	glutAddSubMenu(   "Depth Buffer",  depthbuffermenu);
#endif

#ifdef DEMO_Z_FIGHTING
	glutAddSubMenu(   "Depth Fighting",depthfightingmenu);
#endif

	glutAddSubMenu(   "Depth Cue",     depthcuemenu);
	glutAddSubMenu(   "Projection",    projmenu );
	glutAddMenuEntry( "Reset",         RESET );
	glutAddSubMenu(   "Debug",         debugmenu);
	glutAddMenuEntry( "Quit",          QUIT );

// attach the pop-up menu to the right mouse button:

	glutAttachMenu( GLUT_RIGHT_BUTTON );
}


// the keyboard callback:

void
Keyboard( unsigned char c, int x, int y )
{
	if( DebugOn != 0 )
		fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );

	switch( c )
	{
		case 'o':
		case 'O':
			NowProjection = ORTHO;
			break;

		case 'p':
		case 'P':
			NowProjection = PERSP;
			break;

		case 'q':
		case 'Q':
		case ESCAPE:
			DoMainMenu( QUIT );	// will not return here
			break;				// happy compiler

		default:
			fprintf( stderr, "Don't know what to do with keyboard hit: '%c' (0x%0x)\n", c, c );
	}

	// force a call to Display( ):

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// called when the mouse button transitions down or up:

void
MouseButton( int button, int state, int x, int y )
{
	int b = 0;			// LEFT, MIDDLE, or RIGHT

	if( DebugOn != 0 )
		fprintf( stderr, "MouseButton: %d, %d, %d, %d\n", button, state, x, y );

	
	// get the proper button bit mask:

	switch( button )
	{
		case GLUT_LEFT_BUTTON:
			b = LEFT;		break;

		case GLUT_MIDDLE_BUTTON:
			b = MIDDLE;		break;

		case GLUT_RIGHT_BUTTON:
			b = RIGHT;		break;

		case SCROLL_WHEEL_UP:
			Scale += SCLFACT * SCROLL_WHEEL_CLICK_FACTOR;
			// keep object from turning inside-out or disappearing:
			if (Scale < MINSCALE)
				Scale = MINSCALE;
			break;

		case SCROLL_WHEEL_DOWN:
			Scale -= SCLFACT * SCROLL_WHEEL_CLICK_FACTOR;
			// keep object from turning inside-out or disappearing:
			if (Scale < MINSCALE)
				Scale = MINSCALE;
			break;

		default:
			b = 0;
			fprintf( stderr, "Unknown mouse button: %d\n", button );
	}

	// button down sets the bit, up clears the bit:

	if( state == GLUT_DOWN )
	{
		Xmouse = x;
		Ymouse = y;
		ActiveButton |= b;		// set the proper bit
	}
	else
	{
		ActiveButton &= ~b;		// clear the proper bit
	}

	glutSetWindow(MainWindow);
	glutPostRedisplay();

}


// called when the mouse moves while a button is down:

void
MouseMotion( int x, int y )
{
	int dx = x - Xmouse;		// change in mouse coords
	int dy = y - Ymouse;

	if( ( ActiveButton & LEFT ) != 0 )
	{
		Xrot += ( ANGFACT*dy );
		Yrot += ( ANGFACT*dx );
	}

	if( ( ActiveButton & MIDDLE ) != 0 )
	{
		Scale += SCLFACT * (float) ( dx - dy );

		// keep object from turning inside-out or disappearing:

		if( Scale < MINSCALE )
			Scale = MINSCALE;
	}

	Xmouse = x;			// new current position
	Ymouse = y;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// reset the transformations and the colors:
// this only sets the global variables --
// the glut main loop is responsible for redrawing the scene

void
Reset( )
{
	ActiveButton = 0;
	AxesOn = 1;
	DebugOn = 0;
	DepthBufferOn = 1;
	DepthFightingOn = 0;
	DepthCueOn = 0;
	Scale  = 1.0;
	ShadowsOn = 0;
	NowColor = YELLOW;
	NowProjection = ORTHO;
	Xrot = Yrot = 0.;
}


// called when user resizes the window:

void
Resize( int width, int height )
{
	// don't really need to do anything since window size is
	// checked each time in Display( ):

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// handle a change to the window's visibility:

void
Visibility ( int state )
{
	if( DebugOn != 0 )
		fprintf( stderr, "Visibility: %d\n", state );

	if( state == GLUT_VISIBLE )
	{
		glutSetWindow( MainWindow );
		glutPostRedisplay( );
	}
	else
	{
		// could optimize by keeping track of the fact
		// that the window is not visible and avoid
		// animating or redrawing it ...
	}
}



///////////////////////////////////////   HANDY UTILITIES:  //////////////////////////


// the stroke characters 'X' 'Y' 'Z' :

static float xx[ ] = { 0.f, 1.f, 0.f, 1.f };

static float xy[ ] = { -.5f, .5f, .5f, -.5f };

static int xorder[ ] = { 1, 2, -3, 4 };

static float yx[ ] = { 0.f, 0.f, -.5f, .5f };

static float yy[ ] = { 0.f, .6f, 1.f, 1.f };

static int yorder[ ] = { 1, 2, 3, -2, 4 };

static float zx[ ] = { 1.f, 0.f, 1.f, 0.f, .25f, .75f };

static float zy[ ] = { .5f, .5f, -.5f, -.5f, 0.f, 0.f };

static int zorder[ ] = { 1, 2, 3, 4, -5, 6 };

// fraction of the length to use as height of the characters:
const float LENFRAC = 0.10f;

// fraction of length to use as start location of the characters:
const float BASEFRAC = 1.10f;

//	Draw a set of 3D axes:
//	(length is the axis length in world coordinates)

void
Axes( float length )
{
	glBegin( GL_LINE_STRIP );
		glVertex3f( length, 0., 0. );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., length, 0. );
	glEnd( );
	glBegin( GL_LINE_STRIP );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., 0., length );
	glEnd( );

	float fact = LENFRAC * length;
	float base = BASEFRAC * length;

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 4; i++ )
		{
			int j = xorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( base + fact*xx[j], fact*xy[j], 0.0 );
		}
	glEnd( );

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 5; i++ )
		{
			int j = yorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( fact*yx[j], base + fact*yy[j], 0.0 );
		}
	glEnd( );

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 6; i++ )
		{
			int j = zorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( 0.0, fact*zy[j], base + fact*zx[j] );
		}
	glEnd( );

}


// function to convert HSV to RGB
// 0.  <=  s, v, r, g, b  <=  1.
// 0.  <= h  <=  360.
// when this returns, call:
//		glColor3fv( rgb );

void
HsvRgb( float hsv[3], float rgb[3] )
{
	// guarantee valid input:

	float h = hsv[0] / 60.f;
	while( h >= 6. )	h -= 6.;
	while( h <  0. ) 	h += 6.;

	float s = hsv[1];
	if( s < 0. )
		s = 0.;
	if( s > 1. )
		s = 1.;

	float v = hsv[2];
	if( v < 0. )
		v = 0.;
	if( v > 1. )
		v = 1.;

	// if sat==0, then is a gray:

	if( s == 0.0 )
	{
		rgb[0] = rgb[1] = rgb[2] = v;
		return;
	}

	// get an rgb from the hue itself:
	
	float i = (float)floor( h );
	float f = h - i;
	float p = v * ( 1.f - s );
	float q = v * ( 1.f - s*f );
	float t = v * ( 1.f - ( s * (1.f-f) ) );

	float r=0., g=0., b=0.;			// red, green, blue
	switch( (int) i )
	{
		case 0:
			r = v;	g = t;	b = p;
			break;
	
		case 1:
			r = q;	g = v;	b = p;
			break;
	
		case 2:
			r = p;	g = v;	b = t;
			break;
	
		case 3:
			r = p;	g = q;	b = v;
			break;
	
		case 4:
			r = t;	g = p;	b = v;
			break;
	
		case 5:
			r = v;	g = p;	b = q;
			break;
	}


	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
}

void
Cross(float v1[3], float v2[3], float vout[3])
{
	float tmp[3];
	tmp[0] = v1[1] * v2[2] - v2[1] * v1[2];
	tmp[1] = v2[0] * v1[2] - v1[0] * v2[2];
	tmp[2] = v1[0] * v2[1] - v2[0] * v1[1];
	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}

float
Dot(float v1[3], float v2[3])
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


float
Unit(float vin[3], float vout[3])
{
	float dist = vin[0] * vin[0] + vin[1] * vin[1] + vin[2] * vin[2];
	if (dist > 0.0)
	{
		dist = sqrtf(dist);
		vout[0] = vin[0] / dist;
		vout[1] = vin[1] / dist;
		vout[2] = vin[2] / dist;
	}
	else
	{
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}
	return dist;
}


float
Unit( float v[3] )
{
	float dist = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	if (dist > 0.0)
	{
		dist = sqrtf(dist);
		v[0] /= dist;
		v[1] /= dist;
		v[2] /= dist;
	}
	return dist;
}


float solveDivergence(float p0, float p1, float p2, float p3, float sourceTerm) {
	// Solve Poisson for divergence at each point 
	return 0.25 * (p0 + p1 + p2 + p3 - sourceTerm * CELLSIZE * CELLSIZE);
}

float cubicInterpolate(float p0, float p1, float p2, float p3, float dt) {
	// Perform Centripetal Catmull-Rom interpolation to determine advected value
	// Calculate spline values, use coefficients from Catmull-Rom spline matrix
	float a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
	float a1 = p0 - 2.5 * p1 + 2. * p2 - 0.5 * p3;
	float a2 = -0.5 * p0 + 0.5 * p2;
	// Return the cubic polynomial
	return dt * (dt * (a0 * dt + a1) + a2) + p1;
}

void calculateAdvection(float timestep) {
	// Use the MacCormack Method
	// Do a forward pass to estimate advection at each point using the Semi-Lagrangian method
	//   Trace forward along the velocity field to the estimated new location of the fluid
	//   Interpolate at that point
	//   Store as forward estimate
	// Do a backward pass using Semi-Lagrangian
	//   Trace backward from the new point
	//   Interpolate value at that point
	//   Store as backward estimate
	// Corrected Value = Forward Estimate + 0.5 * (Original Value - Backward Estimate)
	// Then, project onto a divergence-free field using Poisson projection and a Jacobi solver
	// Save the corrected point in the temporary struct, then swap out when done
	// Added in minmod limiter, maybe advance to a MC limiter or Van Leer
	// Currently doing bicubic interpolation for the 4x4 grid [i - 1, j - 1] to [i + 2, j + 2]
	// Bilinear would be more performant but less accurate, add as an option later

	// Create struct of temporary values
	struct Point* tempGrid = (struct Point*)malloc((SIZE * SIZE) * sizeof(Point));

	// Only advect interior cells
	for (int x = 1; x < SIZE - 1; x++) {
		//fprintf(stdout, "%d,,", x);
		for (int y = 1; y < SIZE - 1; y++) {
			float forwardX, forwardY;
			float backwardX, backwardY;
			// Calculate the forward estimate
			// Get next (x, y) coords
			float xNext = (float)x + getAtIndex(x, y).vx * timestep;
			float yNext = (float)y + getAtIndex(x, y).vy * timestep;
			// Round to nearest whole index
			int i = floor(xNext);
			int j = floor(yNext);
			// Calculate fractional offsets
			float dx = xNext - i;
			float dy = yNext - j;
			// Interpolate x and y velocities at the new point
			// vx, vy = CI(CI(p0..p3, dx), dy)
			// If point is outside the grid, use ghost cell
			// Add more quantities (temperature, density, etc) here
			// Maybe turn forward into a forwardStruct at some point?

			// Estimate row-by-row x-velocity
			float tempX[4];
			float tempY[4];
			float vx = getAtIndex(x, y).vx;
			float vy = getAtIndex(x, y).vy;
			for (int n = -1; n <= 2; n++) {
				// Pull 4x4 neighbors, use ghost cells if OOBs
				int u = i + n;
				int v0 = j - 1;
				int v1 = j;
				int v2 = j + 1;
				int v3 = j + 2;

				struct Point p0 = getAtIndex(u, v0);
				struct Point p1 = getAtIndex(u, v1);
				struct Point p2 = getAtIndex(u, v2);
				struct Point p3 = getAtIndex(u, v3);

				tempX[n + 1] = cubicInterpolate(
					p0.vx,
					p1.vx,
					p2.vx,
					p3.vx,
					dx
				);
				tempY[n + 1] = cubicInterpolate(
					p0.vy,
					p1.vy,
					p2.vy,
					p3.vy,
					dx
				);
			}

			// Calculate forward estimates
			forwardX = cubicInterpolate(tempX[0], tempX[1], tempX[2], tempX[3], dy);
			forwardY = cubicInterpolate(tempY[0], tempY[1], tempY[2], tempY[3], dy);

			// Back trace the point and interpolate x and y velocities
			// vx, vy = CI(CI(p0..p3, dx), dy)
			// If point is outside the grid, use ghost cell
			// Add more quantities (temperature, density, etc) here
			// Maybe turn backward into a backStruct at some point?

			float xPrev = (float)x - forwardX * timestep;
			float yPrev = (float)y - forwardY * timestep;
			// Round to nearest whole index
			i = floor(xPrev);
			j = floor(yPrev);
			// Calculate fractional offsets
			dx = xPrev - i;
			dy = yPrev - j;

			float tempXPrev[4];
			float tempYPrev[4];

			for (int n = -1; n <= 2; n++) {
				// Pull 4x4 neighbors, use ghost cells if OOBs
				int u = i + n;
				int v0 = j - 1;
				int v1 = j;
				int v2 = j + 1;
				int v3 = j + 2;
				
				struct Point p0 = getAtIndex(u, v0);
				struct Point p1 = getAtIndex(u, v1);
				struct Point p2 = getAtIndex(u, v2);
				struct Point p3 = getAtIndex(u, v3);

				tempXPrev[n + 1] = cubicInterpolate(
					p0.vx,
					p1.vx,
					p2.vx,
					p3.vx,
					dx
				);
				tempYPrev[n + 1] = cubicInterpolate(
					p0.vy,
					p1.vy,
					p2.vy,
					p3.vy,
					dx
				);
			}

			// Calculate backward estimates
			backwardX = cubicInterpolate(tempXPrev[0], tempXPrev[1], tempXPrev[2], tempXPrev[3], dy);
			backwardY = cubicInterpolate(tempYPrev[0], tempYPrev[1], tempYPrev[2], tempYPrev[3], dy);

			// Calculate corrected value
			//  Corrected Value = Forward Estimate + 0.5 * (Original Value - Backward Estimate)
			float correctedX = forwardX + 0.5 * (getAtIndex(x, y).vx - backwardX);
			float correctedY = forwardY + 0.5 * (getAtIndex(x, y).vy - backwardY);

			// Apply minmod to prevent overshooting
			if (MINMOD) {
				correctedX = minmod(correctedX, forwardX - getAtIndex(x, y).vx);
				correctedY = minmod(correctedY, forwardY - getAtIndex(x, y).vy);
			}
			

			// Copy new value into temporary struct
			// Don't write into the boundary
			if (x > 0 && x < SIZE - 1 && y > 0 && y < SIZE - 1) {
				tempGrid[SIZE * x + y].vx = correctedX;
				tempGrid[SIZE * x + y].vy = correctedY;
			}
			
			if (DEBUG)
			{    // Debug stuff
				//float oldX = getAtIndex(x, y, grid).vx;
				//float oldY = getAtIndex(x, y, grid).vy;
				//float oldMag = sqrt(oldX * oldX + oldY * oldY);
				//float newMag = sqrt(correctedX * correctedX + correctedY * correctedY);
				//fprintf(stdout, "%.3f,", newMag);
			}
		}
		if (DEBUG) {
			//fprintf(stdout, "\n")'
		}
	}

	// Solve for divergence (ensure no areas in the fluid are compressed)
	// All points should have initial pressures = 0
	// Loop k times (figure out what an acceptable divergence looks like and modify later)
	// Initialize pressures to 0 in preparation for Poisson-ification
	float* tempPressures = (float*)calloc(SIZE * SIZE, sizeof(float));
	float* tempOldPressures = (float*)calloc(SIZE * SIZE, sizeof(float));

	// Calculate sourceTerm ((x-comp + y-comp) * material derivative * h^2) for each point
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			float p0 = getAtIndex(i + 1, j).vx;
			float p1 = getAtIndex(i - 1, j).vx;
			float p2 = getAtIndex(i, j + 1).vy;
			float p3 = getAtIndex(i, j - 1).vy;
			tempGrid[SIZE * i + j].sourceTerm = (p0 - p1 + p2 - p3) / (2. * (float)CELLSIZE);
		}
	}
		// Jacobi
		float maxError = INFINITY;
		while (maxError > JACOBIS) {
			// Move over old pressures
			maxError = 0.;
			memcpy(tempOldPressures, tempPressures, SIZE * SIZE * sizeof(float));
			for (int i = 1; i < SIZE - 1; i++) {
				for (int j = 1; j < SIZE - 1; j++) {
					tempPressures[SIZE * i + j] = solveDivergence(
						getSafePressure(tempOldPressures, i + 1, j),
						getSafePressure(tempOldPressures, i - 1, j),
						getSafePressure(tempOldPressures, i, j + 1),
						getSafePressure(tempOldPressures, i, j - 1),
						tempGrid[SIZE * i + j].sourceTerm
					);
					float currentError = fabs(tempOldPressures[i * SIZE + j] - tempPressures[i * SIZE + j]);
					if (currentError > maxError) {
						maxError = currentError;
					}
					fprintf(stderr, "maxError = %f\n", maxError);
				}
			}
			// Update ghost cells
			for (int i = 0; i < SIZE; i++) {
				// top/bottom
				tempPressures[i] = tempPressures[SIZE + i];           // top row
				tempPressures[(SIZE-1) * SIZE + i] = tempPressures[(SIZE-2) * SIZE + i]; // bottom row
				// left/right
				tempPressures[i * SIZE] = tempPressures[i * SIZE + 1];           // left col
				tempPressures[i * SIZE + (SIZE-1)] = tempPressures[i * SIZE + (SIZE-2)]; // right col
			}
		}

		// Now correct x- and y-velocities to account for pressure buildup (or lack thereof)
		// vx = (temp vx) - ((dt / rho) * (Pright - Pleft) / 2 * CELLSIZE)
		// vy = (temp vy) - ((dt / rho) * (Pabove - Pbelow) / 2 * CELLSIZE)

		for (int x = 1; x < SIZE - 1; x++) {
			for (int y = 1; y < SIZE - 1; y++) {
				float scale = timestep / grid[SIZE * x + y].density;
				float h = (float)CELLSIZE;
				float xcorrection = scale * ((getSafePressure(tempPressures, x + 1, y) - getSafePressure(tempPressures, x - 1, y)) / (2. * h));
				float ycorrection = scale * ((getSafePressure(tempPressures, x, y + 1) - getSafePressure(tempPressures, x, y - 1)) / (2. * h));
				tempGrid[x * SIZE + y].vx -= xcorrection;
				tempGrid[x * SIZE + y].vy -= ycorrection;
			}
		}

		// Project temporary grid onto permanent grid
		memcpy(grid, tempGrid, SIZE * SIZE * sizeof(Point));
		free(tempGrid);
		free(tempPressures);
		free(tempOldPressures);
		// Correct boundary cells
		InitBoundaries();
}

struct Point getAtIndex(int x, int y) {
	// Retrieve a point at the given index from the grid
	// Grid[x][y] = grid[SIZE * x + y]
	// If pulling from outside the grid, pull a ghost cell
	x = (x < 0) ? 0 : x;
	x = (x > SIZE - 1) ? SIZE - 1 : x;
	y = (y < 0) ? 0 : y;
	y = (y > SIZE - 1) ? SIZE - 1 : y;
	return grid[SIZE * x + y];
}

void setVxAtIndex(int x, int y, float vx) {
	// Set x-velocity at a given index
	grid[SIZE * x + y].vx = vx;
}

void setVyAtIndex(int x, int y, float vy) {
	grid[SIZE * x + y].vy = vy;
}

float getSafePressure(float* pressureGrid, int x, int y) {
	// Retrieve a point at the given index from the grid
	// Grid[x][y] = grid[SIZE * x + y]
	// If pulling from outside the grid, pull a ghost cell
	x = (x < 0) ? 0 : x;
	x = (x > SIZE - 1) ? SIZE - 1 : x;
	y = (y < 0) ? 0 : y;
	y = (y > SIZE - 1) ? SIZE - 1 : y;
	return pressureGrid[SIZE * x + y];
}

void InitGrid() {
	grid = (struct Point*)malloc((SIZE * SIZE) * sizeof(Point));
	float scale = (float)(SIZE) / 2.;
	if (grid) {
		const float Vmax = 1.f;     // maximum velocity magnitude
		const float Rmin = 0.5f;    // minimum radius to avoid singularity at center

		for (int row = 0; row < SIZE; row++) {
			for (int col = 0; col < SIZE; col++) {
				int idx = SIZE * row + col;

				// Boundary cells
				if (row == 0 || row == SIZE - 1 || col == 0 || col == SIZE - 1) {
					grid[idx].temperature = 0.f;
					grid[idx].vx = 0.f;
					grid[idx].vy = 0.f;
					grid[idx].pressure = 0.f;
					grid[idx].oldPressure = 0.f;
					grid[idx].sourceTerm = 0.f;
					grid[idx].density = 1.f;
					continue;
				}

				// Interior cells: rotational field
				float dx = col - scale;
				float dy = row - scale;
				float r = sqrtf(dx * dx + dy * dy);
				if (r < Rmin) r = Rmin;

				float vx = -dy / r;
				float vy = dx / r;

				// Scale velocity to maximum allowed
				float mag = sqrtf(vx * vx + vy * vy);
				if (mag > Vmax) {
					vx = vx * (Vmax / mag);
					vy = vy * (Vmax / mag);
				}

				grid[idx].temperature = 0.f;
				grid[idx].vx = vx;
				grid[idx].vy = vy;
				grid[idx].pressure = 0.f;
				grid[idx].oldPressure = 0.f;
				grid[idx].sourceTerm = 0.f;
				grid[idx].density = 1.f;
			}
		}
		float maxVelocity = 0.f;
		for (int i = 0; i < SIZE * SIZE; i++) {
			float mag = sqrtf(grid[i].vx * grid[i].vx + grid[i].vy * grid[i].vy);
			if (mag > maxVelocity) maxVelocity = mag;
		}
		printf("Max initial velocity: %f\n", maxVelocity);
	}
	else {
		fprintf(stderr, "Error initializing grid\n");
	}
}

void InitBoundaries() {
	// Each boundary cell should have the negated normal component of its nearest interior cell while keeping the tangential velocity
	// Corners are funky and need the diagonal point

	// Do the top edge [0, SIZE - 1]
	// Set the top left and top right corners to 0
	grid[0].vx = 0.;
	grid[0].vy = 0.;
	grid[SIZE - 1].vx = 0.;
	grid[SIZE - 1].vy = 0.;
	for (int y = 1; y < SIZE - 1; y++) {
		grid[y].vx = getAtIndex(1, y).vx;
		grid[y].vy = -getAtIndex(1, y).vy;
	}

	// Do the bottom edge
	// Set the bottom left and bottom right corners
	grid[(SIZE - 1) * SIZE].vx = 0.;
	grid[(SIZE - 1) * SIZE].vy = 0.;
	grid[(SIZE - 1) * SIZE + (SIZE - 1)].vx = 0.;
	grid[(SIZE - 1) * SIZE + (SIZE - 1)].vy = 0.;
	int x = SIZE - 1;
	for (int y = 1; y < SIZE - 1; y++) {
		grid[x * SIZE + y].vx = grid[(x - 1) * SIZE + y].vx;
		grid[x * SIZE + y].vy = -grid[(x - 1) * SIZE + y].vy;
	}

	// Do the left side
	for (int x = 1; x < SIZE - 1; x++) {
		grid[x * SIZE].vx = -grid[x * SIZE + 1].vx;
		grid[x * SIZE].vy = grid[x * SIZE + 1].vy;
	}

	// Do the right side
	int y = SIZE - 2;
	for (int x = 1; x < SIZE - 1; x++) {
		grid[x * SIZE + SIZE - 1].vx = -grid[x * SIZE + y].vx;
		grid[x * SIZE + SIZE - 1].vy = grid[x * SIZE + y].vy;
	}
}

float minmod(float a, float b) {
	// Take more conservative slope if a and b have same sign
	// Flatten to 0 if opposite signs
	if (a * b <= 0.0) {
		return 0.0;
	}
	return (fabs(a) < fabs(b)) ? a : b;
}


uint64_t timeSinceEpochMillisec() {
	using namespace std::chrono;
	return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}
