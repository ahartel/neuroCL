// Pong by Konstantinos Egarhos (aka konsnos)
/**
* ...
* @author konsnos
*/

#include <GL/freeglut.h>	// OpenGL toolkit
#include <stdio.h>
#include "compare/network.h"

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
GLdouble bpos_old;
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

/*********************************************
    Neuron connections
*********************************************/
void updateNeurons()
{
	std::vector<pre_spike> spikes;

//	for (size_t ii=0;ii<Height; ++ii)
//		neurons[ii] = 0;

	if (p2_old > p2)
		for (size_t ii=p2/2;ii<p2_old/2; ++ii)
			spikes.push_back(pre_spike(ii,0,10));
			//neurons[ii] = 1;
	else
		for (size_t ii=p2_old/2;ii<p2/2; ++ii)
			spikes.push_back(pre_spike(ii,0,10));
			//neurons[ii] = 1;

	if (bpos_old > bpos[1])
		for (size_t ii=bpos[1]/2;ii<bpos_old/2; ++ii)
			spikes.push_back(pre_spike(ii+250,0,10));
			//neurons[ii+250] = 1;
	else
		for (size_t ii=bpos_old/2;ii<bpos[1]/2; ++ii)
			spikes.push_back(pre_spike(ii+250,0,10));
			//neurons[ii+250] = 1;

	//spikes
	//for (size_t ii=0;ii<Height; ++ii)
	//{

	network.add_spikes(spikes);
	network.step();
	int direction = 0;
	std::vector<unsigned int> get_spikes = network.get_last_spikes();

	if (!get_spikes.empty())
	{
		for (auto sp : get_spikes)
		{
			if (sp>=790 && sp<800)
				direction += 1;
			else if (sp>=780 && sp<790)
				direction -= 1;
		}
	}

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
	 bpos_old = bpos[1];
     // Ball collision with the border of the window.
     if(bpos[0] + bspeed > 500-bsize)
     {
          bvx = -1;
          score[0]++;
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
     if(bpos[0] <= 10 + paddleWidth && bpos[1] >= p1 && bpos[1] + bsize <= p1 + paddleHeight)
     {
          bvx *= -1;
          // bspeed *= 1.01;
     }
     // Collision with right paddle.
     if(bpos[0] + bsize >= 480 && bpos[1] >= p2 && bpos[1] + bsize <= p2 + paddleHeight)
     {
          bvx *= -1;
          // bspeed *= 1.01;
     }
     
     // New position of the ball dependent to the speed ***************
     // Because frame rate isn't always the same that might affect objects speed.
     // That's why we multiply the frame speed and multiply it with the object speed.
     // In that way the object speed is always the same.
     bpos[0] += bspeed * bvx * runtime;
     bpos[1] += bspeed * bvy * runtime;

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

	 updateNeurons();
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
	if (curTime < 10000)
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
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

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

