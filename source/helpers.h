#pragma once
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>

using namespace std;

/*
 * Constants
 *
 * Good values for debugging are Ne=4, Ni=1, M=6
 */
const int Ne = 800;//12800;
const int Ni = 200;//3200;
const int N = Ne+Ni; // total number of neurons
const static int M = 100; // number of postsynaptic neurons
const float v_thresh = 30.0/*mV*/;
const float v_reset = -65.0/*mV*/;
const int T = 100; // number of seconds to simulate
const float h = 1; // timestep in milliseconds
const static int D = 20; // max. delay in ms
const float	sm = 10.0;		// maximal synaptic strength

const unsigned int neurons_tobe_watched[] = {0,1,2,3,4};
const unsigned int num_watched_neurons = 5;


#define getrandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1

#define TIMING
#undef DEBUG_OUTPUT
#undef WATCH_NEURONS
#undef WATCH_DERIVATIVES
#undef WATCH_ADAPTATION

template<typename T>
void print_loop (T array, size_t size)
{
	if (size > 0)
	{
		if (array[0] < 10)
			cout << 0;
		cout << array[0];
		for (unsigned int i=1; i<size; i++)
		{
			cout << ", ";
			if (array[i] < 10)
				cout << 0;
			cout << array[i];
		}
	}
	cout << endl;
}

template<typename T>
void print_loop (T array)
{
	print_loop(array,array.size());
}

void write_derivatives(string filename, unsigned int sec, std::vector<std::vector<float> >const& watched_LTP, std::vector<std::vector<float> > const& watched_LTD);

void write_spikes(string filename, unsigned int sec, vector<unsigned int> const& spikes, unsigned int* k);

void write_watched_membranes(string filename, unsigned int sec, std::vector<std::vector<float> >const& watched_membrane, std::vector<std::vector<float> > const& watched_us);

void init_neurons_sparse(
	vector<float> & membranes,
	vector<float> & u,
	vector<float> & d,
	vector<float> & a,
	float I[N][D+1],
	vector<float> weights[N],
	unsigned int delay_start[N][D],
	unsigned int delay_count[N][D],
	vector<unsigned int> post_neurons[N],
	unsigned int num_post[N],
	float LTP[N][1001+D],
	float LTD[N],
	vector< vector<float> >& sd,
	vector< vector<float*> >& sd_pre,
	unsigned int num_pre[N],
	vector<unsigned int> pre_neurons[N],
	unsigned int delay_start_pre[N][D],
	unsigned int delay_count_pre[N][D]
	);

void init_neurons(float* membranes, float* u, float* d, float* a, float I[N][D], float weights[N][M], int delays[N][M], int post_neurons[N][M], int* num_post);

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



