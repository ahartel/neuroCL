#include <stdlib.h>
#include <iostream>
#include <vector>
/*
 * Constants
 */
const int Ne = 800;//12800;
const int Ni = 200;//3200;
const int N = Ne+Ni; // total number of neurons
const static int M = 100; // number of postsynaptic neurons
const float v_thresh = 30./*mV*/;
const float v_reset = -65/*mV*/;
const int T = 100; // number of seconds to simulate
const float h = 1; // timestep in milliseconds
const static int D = 20; // max. delay in ms

using namespace std;

#define getrandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1

#define TIMING
#undef WATCH_NEURONS
#undef WATCH_ADAPTATION

void init_neurons_sparse(float* membranes, float* u, float* d, float* a, float I[N][D], vector<float> weights[N], unsigned int delay_start[N][D], unsigned int delay_count[N][D], vector<unsigned int> post_neurons[N], unsigned int num_post[N])
{
	srand(42);

	/*
	 * init neuron parameters
	 */
	// iterate over neurons
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
	}

	/*
	 * init connectivity
	 */
	// iterate over neurons
	for (int i=0; i<N; i++) {
		num_post[i] = getrandom(M);
		int delays[M];
		vector<int> synapse_numbers[D];
		vector<int> pst_n;
		vector<float> wgt_n;
		// set delay_count and delay_start initially to zero
		for (int d=0; d<D; d++)
		{
			delay_start[i][d] = 0;
			delay_count[i][d] = 0;
		}
		// create values for connectivity
		for (int m=0; m<num_post[i]; m++)
		{
			if (i<Ne)
				wgt_n.push_back(6.);
			else
				wgt_n.push_back(-5.);

			delays[m] = getrandom(D);
			synapse_numbers[delays[m]].push_back(m);
			delay_count[i][delays[m]]++;

			pst_n.push_back(getrandom(N));

			while (pst_n.back() == i)
				pst_n[pst_n.size()-1] = getrandom(N);

			//cout << "Neuron " << i << ": Added post-neuron " << pst_n.back() << " with delay " << delays[m] << endl;
		}
		// sort connectivity values according to delay values
		vector<unsigned int> pst_n_sort;
		vector<float> wgt_n_sort;
		for (int d=0; d<D; d++)
		{
			vector<int>::iterator delay_it = synapse_numbers[d].begin();
			while (delay_it != synapse_numbers[d].end())
			{
				pst_n_sort.push_back(pst_n[*delay_it]);
				wgt_n_sort.push_back(wgt_n[*delay_it]);

				delay_it++;
			}
		}

		// add values to global arrays
		weights[i] = wgt_n_sort;
		post_neurons[i] = pst_n_sort;

		int last_delay = 0;
		for (int d=0; d<D; d++)
		{
			if (delay_count[i][d] > 0)
				delay_start[i][d] = last_delay;
			else
				delay_start[i][d] = 0;
			last_delay += delay_count[i][d];
		}
	}
}

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

		num_post[i] = getrandom(M);

		for (int m=0; m<num_post[i]; m++)
		{
			if (i<Ne)
				weights[i][m] = 6.;
			else
				weights[i][m] = -5.;

			delays[i][m] = getrandom(D);

			post_neurons[i][m] = getrandom(N);
			while (post_neurons[i][m] == i)
				post_neurons[i][m] = getrandom(N);
		}
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

const unsigned int neurons_tobe_watched[] = {0,17};
const unsigned int num_watched_neurons = 2;

