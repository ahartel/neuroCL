/*
 * Constants
 */
const int Ne = 800;//12800;
const int Ni = 200;//3200;
const int N = Ne+Ni; // total number of neurons
const int M = 1; // number of postsynaptic neurons
const float v_thresh = 30./*mV*/;
const float v_reset = -65/*mV*/;
const int T = 100000; // number of timesteps
const float h = 0.001; // timestep in seconds

void init_neurons(float* membranes, float* u, float* d, float* a, float* I)
{
	for (int i=0; i<N; i++) {
		//membranes[i] = (float)rand()/float(RAND_MAX)*50.0;
		membranes[i] = v_reset;
		//I[i] = 10.0;
		I[i] = (float)rand()/float(RAND_MAX)*150.0;
		u[i] = 0.2*membranes[i];
		if (i < Ne) {
			a[i] = 0.02;
			d[i] = 8.0;
		}
		else {
			a[i] = 0.1;
			d[i] = 2.0;
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
