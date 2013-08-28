// Pong by Konstantinos Egarhos (aka konsnos)
/**
* ...
* @author konsnos
*/

#include <GL/freeglut.h>	// OpenGL toolkit
#include <stdio.h>
#include "compare/network.h"
#include <fstream>

using namespace std;

// Function Initialising
void drawRect(GLfloat, GLfloat, int, int, float, float, float);
void update();

// Window dimensions
const int Width = 500;
const int Height = 500;

// Neurons
Network network;
//bool neurons[Height];

// Dimensions
#define paddleHeight 80
#define paddleWidth 10
#define zAxis 0

// lower y-values of paddles
float p1 = 220;
float p2 = 220;
float p2_old;
// paddle speeds
float p1speed = 0;
float p1speedup = 3;
float p2speed = 0;
float p2speedup = 3;
// GLfloat player1[4][2] = { 10, p1, 10, p1+paddleHeight, 10+paddleWidth, p1+paddleHeight, 10+paddleWidth, p1 };
// GLfloat player2[4][2] = { 480, p2, 480, p2+paddleHeight, 480+paddleWidth, p2+paddleHeight, 480+paddleWidth, p2 };
const float pspeed = 7;

int score[2];// = {0, 0}; // Score for players 1 and 2

#define bsize 10
GLdouble bpos_old[2] = {235,235};
GLdouble bpos[2] = {235,235};
GLfloat bvx = 1; // ball vector in left-right axis
GLfloat bvy = 1; // ball vector in top-bottom axis
// GLfloat ball[4][2] = { pos[0], pos[1], pos[0], pos[1]+bsize, pos[0]+bsize, pos[1]+bsize, pos[0]+bsize, pos[1] };
float bspeed = 3;

// Time variables for the speed calculation of the objects
GLdouble curTime;
GLdouble prtime;
GLdouble timeDiff;
GLdouble runtime;

// Text variable
static char str[200];
void * font = GLUT_BITMAP_9_BY_15;

// performance measurement
unsigned int good = 0;
unsigned int baad = 0;

/*********************************************
    Neuron connections
*********************************************/
void updateNeurons(bool collision)
{
	std::vector<pre_spike> spikes;

//	for (size_t ii=0;ii<Height; ++ii)
//		neurons[ii] = 0;
/*
	if (p2_old > p2)
		for (size_t ii=p2/2;ii<p2_old/2; ++ii)
			spikes.push_back(pre_spike(ii+500,0,10));
	else
		for (size_t ii=p2_old/2;ii<p2/2; ++ii)
			spikes.push_back(pre_spike(ii,0,10));

	if (bpos_old[0] > bpos[0])
		for (size_t ii=bpos[0]/2;ii<bpos_old[0]/2; ++ii)
			spikes.push_back(pre_spike(ii,0,10));
	else
		for (size_t ii=bpos_old[0]/2;ii<bpos[0]/2; ++ii)
			spikes.push_back(pre_spike(ii,0,10));

	if (bpos_old[1] > bpos[1])
		for (size_t ii=bpos[1]/2;ii<bpos_old[1]/2; ++ii)
			spikes.push_back(pre_spike(ii+250,0,10));
	else
		for (size_t ii=bpos_old[1]/2;ii<bpos[1]/2; ++ii)
			spikes.push_back(pre_spike(ii+250,0,10));
*/

	spikes.push_back(pre_spike(bpos[0]/2,0,10));
	spikes.push_back(pre_spike(bpos[1]/2+250,0,10));
	spikes.push_back(pre_spike(p2/2+500,0,10));


	unsigned int sec_790 =0;
	unsigned int sec_791 =0;
	unsigned int fir_790 =0;
	unsigned int fir_791 =0;

	network.add_spikes(spikes);
	if (collision)
		network.eject_dopamine();
	network.step();
	std::vector<unsigned int> first_spikes = network.get_last_spikes();
	network.step();
	std::vector<unsigned int> next_spikes = network.get_last_spikes();
	first_spikes.insert(first_spikes.end(),next_spikes.begin(),next_spikes.end());


	network.step();
	next_spikes = network.get_last_spikes();
	first_spikes.insert(first_spikes.end(),next_spikes.begin(),next_spikes.end());
	network.step();
	next_spikes = network.get_last_spikes();
	first_spikes.insert(first_spikes.end(),next_spikes.begin(),next_spikes.end());

	if (!first_spikes.empty())
	{
		for (auto sp : first_spikes)
		{
			if (sp>=770 && sp<772)
				fir_790 += 1;
			else if (sp>=780 && sp<782)
				fir_791 += 1;
		}
	}

	//if (!first_spikes.empty())
	//{
	//	for (auto sp : first_spikes)
	//	{
	//		if (sp>=770 && sp<785)
	//			sec_790=1;
	//		if (sp>=785 && sp<800)
	//			sec_791=1;
	//	}
	//}


	int direction = 0;

	if (fir_790 > fir_791)
		direction = 1;
	else if (fir_791 > fir_790)
		direction = -1;


	//cout << "up: " << up_spikes << "\tdown: " << down_spikes << endl;

	/*
	bool flip = rand()%2;
	if (flip)
	{
		bool up = rand()%2;
		if (up)
			direction += 1;
		else
			direction += 1;
	}
	*/

	if (direction > 0)
	{
		p2speed = pspeed*p2speedup;
		p2speedup *= 0.9;
	}
	else if (direction < 0)
	{
		p2speed = -pspeed*p2speedup;
		p2speedup *= 0.9;
	}

}

void write_performance()
{
	cout << "Writing performance" << endl;
	ofstream myfile( "performance.txt", ios::app);
	myfile << float(good)/float(baad) << "\n";
	myfile.close();

	good = 0;
	baad = 0;
}

/*********************************************
    Drawing
*********************************************/
void renderScene(void)
{
     // Elapsed time from the initiation of the game.
     curTime = glutGet(GLUT_ELAPSED_TIME);
     timeDiff = curTime - prtime; // Elapsed time from the previous frame.
     prtime = curTime;
     runtime = timeDiff / 10.0;
     
     update(); // Update the game.
     // Clear the screen.
     glClear(GL_COLOR_BUFFER_BIT);
     
     
     sprintf(str, "Player 1 Score: %d, Player 2 Score: %d", score[0], score[1] );
     glRasterPos2f(10, 10);
     glutBitmapString(font,(unsigned char*)str);
    
     // Draws the 1st paddle and the 2nd.
     drawRect(10.0, p1, paddleHeight, paddleWidth, 0, 0, 1); // Left blue paddle.
     drawRect(480.0, p2, paddleHeight, paddleWidth, 1, 0, 0); // Right red paddle.
	
	 drawRect(bpos[0], bpos[1], bsize, bsize, 1, 1, 0); // Yellow ball.
	
     glutSwapBuffers(); // Draw the new frame of the game.
}

/*********************************************
    Update
*********************************************/
void update()
{
	bool right_collision = 0;

     bpos_old[0] = bpos[0];
     bpos_old[1] = bpos[1];
     // Ball collision with the border of the window.
     if(bpos[0] + bspeed > 500-bsize)
     {
          bvx = -1;
          score[0]++;
		  baad += 1;
     }
     else if (bpos[0] - bspeed < 0)
     {
          bvx = 1;
          score[1]++;
     }
     if(bpos[1] + bspeed > 500-bsize)
          bvy = -1;
     else if (bpos[1] - bspeed < 0)
          bvy = 1;

     /******** Collision with the paddles ********/
     // Collision with left paddle.
     if(bpos[0] <= 10 + paddleWidth && bpos[1] >= p1 && bpos[1] + bsize <= p1 + paddleHeight && bvx == -1)
     {
          bvx *= -1;
          // bspeed *= 1.01;
     }
     // Collision with right paddle.
     if(bpos[0] + bsize >= 480 && bpos[1] >= p2 && bpos[1] + bsize <= p2 + paddleHeight && bvx == 1)
     {
          bvx *= -1;
          right_collision = 1;
		  good += 1;
          // bspeed *= 1.01;
     }
     
     // New position of the ball dependent to the speed ***************
     // Because frame rate isn't always the same that might affect objects speed.
     // That's why we multiply the frame speed and multiply it with the object speed.
     // In that way the object speed is always the same.
     bpos[0] += bspeed * bvx;// * runtime;
     bpos[1] += bspeed * bvy;// * runtime;

	 // Player 1
	 if (p1speed > 0.1 || p1speed < -0.1)
	 {
		 if(420>=p1+p1speed && 0<=p1+p1speed)
		 {
			 p1 += p1speed;
			 p1speed = p1speed * 0.8;
		 }
	 }
	 else
		 p1speed = 0;
	 if (p1speedup <= 1.8)
		 p1speedup += 0.2;
	 // Player 2
	 p2_old = p2;
	 if (p2speed > 0.1 || p2speed < -0.1)
	 {
		 if(420>=p2+p2speed && 0<=p2+p2speed)
		 {
			 p2 += p2speed;
			 p2speed = p2speed * 0.8;
		 }
	 }
	 else
		 p2speed = 0;
	 if (p2speedup <= 1.8)
		 p2speedup += 0.2;

	 // killer KI
	 p1 = bpos[1];

	 updateNeurons(right_collision);

	 if (good + baad == 20)
		 write_performance();
}


/*********************************************
    Create and draw a rectangle.
    Generally, we are only interested in the height(y) of the paddle.
*********************************************/
void drawRect(GLfloat x, GLfloat y, int height, int width, float R, float G, float B)
{
     glColor3f(R,G,B);
     
     glBegin(GL_QUADS);
          glVertex2f( x, y);				        // Top Left
          glVertex2f( x, y + height);				// Bottom Left
          glVertex2f( x + width, y + height);		// Bottom Right
          glVertex2f( x + width, y);				// Top Right
	 glEnd();
}


/*********************************************
    User Input (keyboard)
    
    I don't calculate frame speed, since it doesn't affect much.
*********************************************/
void keyPress(unsigned char key, int x, int y)
{
    switch (key) 
    {
        case 'w':
        case 'W':
				p1speed = pspeed*p1speedup;
				p1speedup *= 0.9;

             break;
        case 's':
        case 'S':
				p1speed = -pspeed*p1speedup;
				p1speedup *= 0.9;
             break;
        case 27 :
             exit(0);
    }
}

/*********************************************
    Special Keys Input (keyboard arrows)
    
    I don't calculate frame speed, since it doesn't affect much.
*********************************************/
void specialKeyPress(int key, int x, int y)
{
	/*
     switch (key) 
     {
            case GLUT_KEY_UP:
				p2speed = pspeed*p2speedup;
				p2speedup *= 0.9;
                 break;
            case GLUT_KEY_DOWN:
				p2speed = -pspeed*p2speedup;
				p2speedup *= 0.9;
                 break;
     }
*/
}

/*********************************************
    Let there be initial variables.
*********************************************/
void init(void)
{
     // The color the windows will redraw. Its done to erase the previous frame.
     glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black, no opacity(alpha).
     
     // Orthographic projection matrix (initiating 2D window).
     gluOrtho2D(0,500,0,500);
     
     // Initialize time.
     prtime = curTime = glutGet(GLUT_ELAPSED_TIME);   
     
     // Initialize score.
     score[0] = score[1] = 0;
}

/* Main function */
int main(int argc, char *argv[])
{
    // Initialize openGL with Double buffer and RGB color without transparency.
    // Its interesting to try GLUT_SINGLE instead of GLUT_DOUBLE.
    glutInit(&argc, argv);
    //glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitDisplayMode( GLUT_RGB);

    // Create the window.
    glutInitWindowSize(Width, Height);
    glutInitWindowPosition(150,50);
    glutCreateWindow("Pong by kon_nos");

    // Define the draw function.
    glutDisplayFunc(renderScene);

    // Define the keyboard input function.
    glutKeyboardFunc(keyPress);    
    // Define the keyboards' special keys input function.
    glutSpecialFunc(specialKeyPress);
    
    // Define which function to execute when nothing except the gameplay is updating.
    glutIdleFunc(renderScene);

    init();
    
    glutMainLoop();

    return 0;
}

