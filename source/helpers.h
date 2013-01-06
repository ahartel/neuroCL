#include <stdlib.h>
#include <iostream>
/*
 * Constants
 */
const int Ne = 800;//12800;
const int Ni = 200;//3200;
const int N = Ne+Ni; // total number of neurons
const static int M = 100; // number of postsynaptic neurons
const float v_thresh = 30./*mV*/;
const float v_reset = -65/*mV*/;
const int T = 10000; // number of timesteps
const float h = 1; // timestep in milliseconds
const static int D = 20; // max. delay in ms

using namespace std;

#define getrandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1

#define TIMING
#define WATCH_NEURONS
#undef WATCH_ADAPTATION

void init_neurons(float* membranes, float* u, float* d, float* a, float I[N][D], float weights[N][M], int delays[N][M], int post_neurons[N][M], int* num_post)
{
	srand(42);
	for (int i=0; i<N; i++) {
		//membranes[i] = (float)rand()/float(RAND_MAX)*50.0;
		membranes[i] = v_reset;

		for (int d=0; d<D; d++)
		{
			I[i][d] = 0;//(float)rand()/float(RAND_MAX)*10.0;
		}
		u[i] = 0.2*membranes[i];
		if (i < Ne) {
			a[i] = 0.02;
			d[i] = 8.0;
		}
		else {
			a[i] = 0.1;
			d[i] = 2.0;
		}

		num_post[i] = getrandom(100);

		for (int m=0; m<num_post[i]; m++)
		{
			if (i<Ne)
				weights[i][m] = 6.;
			else
				weights[i][m] = -5.;
			delays[i][m] = 1;
			post_neurons[i][m] = getrandom(N);
			while (post_neurons[i][m] == i)
				post_neurons[i][m] = getrandom(N);
		}
		/*
		exc_input[i] = 0.0;

		for (int j=0; j<N; j++) {
			int r;
			do{
				exists = 0;		// avoid multiple synapses
				if (i<Ne) r = getrandom(N);
				else	  r = getrandom(Ne);// inh -> exc only
				if (r==i) exists=1;									// no self-synapses 
				for (int k=0;k<j;k++) if (post[i][k]==r) exists = 1;	// synapse already exists  
			}while (exists == 1);
			post[i][j]=r;
		}
		*/
	}
}

#include <chrono>
#include <iostream>

#ifdef TIMING
#define INIT_TIMER(var) auto var = std::chrono::high_resolution_clock::now();
#define START_TIMER(var)  var = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name,var)  std::cout << "RUNTIME of " << name << ": " << \
std::chrono::duration_cast<std::chrono::microseconds>( \
		std::chrono::high_resolution_clock::now()-var \
).count() << " us " << std::endl;
#else
#define INIT_TIMER(var)
#define START_TIMER(var)
#define STOP_TIMER(name,var)
#endif

#ifdef WATCH_NEURONS
const unsigned int neurons_tobe_watched[] = {0,17};
const unsigned int num_watched_neurons = 2;
#endif

